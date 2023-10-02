#include "AlevinHash.hpp"
#include "AlevinTypes.hpp"
#if 0
#include <unordered_set>
#include <unordered_map>
//#include "spdlog/spdlog.h"
//#include "AlevinOpts.hpp"
#include "AlevinUtils.hpp"
//#include "SingleCellProtocols.hpp"
#include "GZipWriter.hpp"
#include "TranscriptGroup.hpp"
//#include "CollapsedCellOptimizer.hpp"
#include "BarcodeGroup.hpp"
//#include "SalmonIndex.hpp"
#endif

size_t readBfh(bfs::path& eqFilePath,
               std::vector<std::string>& txpNames,
               size_t bcLength,
               EqMapT &countMap,
               std::vector<std::string>& bcNames,
               CFreqMapT& freqCounter,
               TrueBcsT& trueBarcodes,
               size_t& totalNormalized,
               bool hasWhitelist
               ) {
  if (hasWhitelist && trueBarcodes.size() == 0) {
    fmt::print(stderr, "whitelist file had 0 CB");
    std::cerr<<std::endl<<std::flush;;
    return 0;
  }

  std::ifstream equivFile(eqFilePath.string());

  size_t numReads{0};
  size_t numTxps, numBcs, numEqclasses;

  // Number of transcripts
  equivFile >> numTxps;

  // Number of barcodes
  equivFile >> numBcs;

  // Number of equivalence classes
  equivFile >> numEqclasses;

  txpNames.resize(numTxps);
  for (size_t i=0; i<numTxps; i++) {
    equivFile >> txpNames[i] ;
  }

  bcNames.resize(numBcs);
  for (size_t i=0; i<numBcs; i++) {
    equivFile >> bcNames[i] ;

    if (bcNames[i].size() != bcLength) {
      fmt::print(stderr, "CB {} has wrong length", bcNames[i]);
      std::cerr<<std::endl<<std::flush;;
      return 0;
    }
  }

  spp::sparse_hash_map<size_t, size_t> barcodeMap;
  { // bcname and bc index rearrangement based on external whitelist
    if (not hasWhitelist) {
      for (size_t i=0; i<bcNames.size(); i++) {
        barcodeMap[i] = i;
      }
    } else {
      // convert set to indexed vector
      size_t idx{0};
      spp::sparse_hash_map<std::string, size_t> trueBarcodeMap;
      for (auto& bc: trueBarcodes) {
        trueBarcodeMap[bc] = idx;
        idx += 1;
      }

      // extracting relevant barcodes
      for (size_t i=0; i<bcNames.size(); i++) {
        if (trueBarcodeMap.contains(bcNames[i])) {
          barcodeMap[i] = trueBarcodeMap[ bcNames[i] ];
        }
      }

      bcNames.clear();
      bcNames.resize(trueBarcodeMap.size());
      for(auto& it: trueBarcodeMap) {
        bcNames[it.second] = it.first;
      }

    } // end else case of not hasWhitelist
  } // end name/index rearrangement

  countMap.max_num_worker_threads(1);
  countMap.reserve(1000000);

  alevin::types::AlevinUMIKmer umiObj;
  //printing on screen progress
  const char RESET_COLOR[] = "\x1b[0m";
  char green[] = "\x1b[30m";
  green[3] = '0' + static_cast<char>(fmt::GREEN);
  char red[] = "\x1b[30m";
  red[3] = '0' + static_cast<char>(fmt::RED);
  std::cerr<<std::endl;

  for (size_t i=0; i<numEqclasses; i++) {
    uint32_t count;
    size_t labelSize ;
    equivFile >> labelSize;

    std::vector<uint32_t> txps(labelSize);
    for (auto& tid : txps) { equivFile >> tid; }
    auto txGroup = TranscriptGroup (txps);

    size_t bgroupSize;
    equivFile >> count >> bgroupSize;

    size_t normalizer{0};
    uint32_t countValidator {0};
    for (size_t j=0; j<bgroupSize; j++){
      size_t ugroupSize;
      std::string bcName;
      bool skipCB {false};
      uint32_t old_bc, new_bc;

      equivFile >> old_bc >> ugroupSize;
      if (not barcodeMap.contains(old_bc)) {
        skipCB = true;
      } else {
        new_bc = barcodeMap[old_bc];
        bcName = bcNames[new_bc];
      }

      for (size_t k=0; k<ugroupSize; k++){
        std::string umiSeq;
        uint32_t umiCount;

        equivFile >> umiSeq >> umiCount;
        countValidator += umiCount;
        if (skipCB) {normalizer += umiCount; continue;}

        uint64_t umiIndex;
        bool isUmiIdxOk = umiObj.fromChars(umiSeq);
        if(isUmiIdxOk){
          umiIndex = umiObj.word(0);
          auto upfn = [new_bc, umiIndex, umiCount](SCTGValue& x) -> void {
            // update the count
            x.count += umiCount;
            x.updateBarcodeGroup(new_bc, umiIndex, umiCount);
          };

          SCTGValue value(umiCount, new_bc, umiIndex, true);
          countMap.upsert(txGroup, upfn, value);
          freqCounter[bcName] += umiCount;
        }
      }// end-ugroup for
    }//end-bgroup for

    if (count != countValidator){
      fmt::print(stderr, "\nBFH eqclass count mismatch:\t"
                 "{} Orignial, validator {} "
                 "Eqclass number {}",
                 count, countValidator, i);
      std::cerr<<std::endl<<std::flush;;
      return 0;
    }

    totalNormalized += normalizer;
    numReads += (countValidator - normalizer);
    double completionFrac = i*100.0/numEqclasses;
    uint32_t percentCompletion {static_cast<uint32_t>(completionFrac)};
    if ( percentCompletion % 10 == 0 || percentCompletion > 95) {
      fmt::print(stderr, "\r{}Done Reading : {}{}%{} and skipped reads: {}{}{}",
                 green, red, percentCompletion, green,
                 red, normalizer, RESET_COLOR);
    }
  }//end-eqclass for
  std::cerr<<std::endl;
  equivFile.close();

  return numReads;
}
