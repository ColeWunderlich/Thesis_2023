#include "DedupUMI.hpp"
#include "tsl/hopscotch_map.h" //I don't think this is needed, included at top of DedupUMI.cpp, code seems to run fine without it   


#include <iostream>
//#include <GZipWriter.hpp>
#include "CollapsedCellOptimizer.hpp" //this includeds GZipWriter and DedupUMI headers 
//#include "DedupUMI.hpp"

//Adding readBfh() from AlevinHash.cpp
//#include "AlevinHash.hpp"
//#include "AlevinHash.cpp"  //Neeed to make sure cmake can find this ok..
#include "readBFH_only.hpp"

//for debug 
#include "printMethods.hpp"

#define TAKE_PATHS
#define TXPMAP_OVERRIDE
#define CHECK_UNIQ_TXP
#define OFF
//#define MY_DEBUG

void getTxpToGeneMap(spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                      std::vector<std::string>& transcripts,
                      const std::string& geneMapFile,
  		      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap){

   std::ifstream t2gFile(geneMapFile);

   spp::sparse_hash_map<std::string, uint32_t> txpIdxMap(transcripts.size());

   for (size_t i=0; i<transcripts.size(); i++){
     txpIdxMap[ transcripts[i] ] = i;
   }

   uint32_t tid, gid, geneCount{0};
   std::string tStr, gStr;
   if(t2gFile.is_open()) {
     while( not t2gFile.eof() ) {
       t2gFile >> tStr >> gStr;

       if(not txpIdxMap.contains(tStr)){
         continue;
       }
       tid = txpIdxMap[tStr];

       if (geneIdxMap.contains(gStr)){
         gid = geneIdxMap[gStr];
       }
       else{
         gid = geneCount;
         geneIdxMap[gStr] = gid;
         geneCount++;
       }

       txpToGeneMap[tid] = gid;
     }
     t2gFile.close();
   }
   if(txpToGeneMap.size() < transcripts.size()){
     std::cerr << "ERROR: "
               << "Txp to Gene Map not found for "
               << transcripts.size() - txpToGeneMap.size()
               <<" transcripts. Exiting" << std::flush;
     exit(1);
   }
 }



//multi-threaded portion of umihandling 
void mtDedupUMIPerCB( EqMapT& fullEqMap, 
		      std::vector<std::string>& bcNames,
		      std::atomic<uint32_t>& barcode,
		      std::deque<std::pair<TranscriptGroup, uint32_t>>& orderedTgroups,
		      std::vector<uint32_t>& umiCounts,
		      size_t totalCells,
		      size_t totalGenes,
		      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
		      GZipWriter& gzw,
		      std::shared_ptr<spdlog::logger>& jointlog,
		      uint32_t umiLength,
		      uint32_t umiEditDistance,
		      bool dumpUmiGraph, 
		      bool dumpArborescences,
		      bool verbose,
		      bool quiet,
		      std::atomic<uint64_t>& totalUniEdgesCount,
             	      std::atomic<uint64_t>& totalBiEdgesCount,
		      std::mutex& dedupMtx,
		      std::ofstream& dedupStream)
{
	//size_t numCells {trueBarcodes.size()}; //for progress at end
  	size_t barcodeIdx;


	//NOTE: Consider removing fragment count validator as it is redundant

	//Run thread untill there are no cells left to process
	//++ on int atomic is equivalent to .fetch_add(1) with sequential consistent memory order (strictest)
	//Every loop starts with previous value of barcode (ie (b=a++) is equivalent to a's old value)
	while((barcodeIdx = barcode++) < totalCells) 
	{
    		//skip the barcode if no mapped UMI
    		if ( umiCounts[barcodeIdx] == 0 ) 
		{
    		  //skippedCB[trueBarcodeIdx].inActive = true; //skipped cell count turned off for now 
    		  continue;
    		}

    		//extract CB string 
    		auto& barcodeStr = bcNames[barcodeIdx];
		
		double totalCount{0.0};
        	//double totalExpGenes{0};
        	std::vector<UGroupT> umiGroups;
        	std::vector<tgrouplabelt> txpGroups;
        	std::vector<double> geneAlphas(totalGenes, 0.0);
        	std::vector<uint8_t> tiers (totalGenes, 0);
		
		//Only for verbose=T (CB+Eqc+TotalUMICounts)
        	std::vector<uint32_t> eqIDs;
        	std::vector<uint32_t> eqCounts;
	
		//process
		size_t fragmentCountValidator {0};
		for (auto& key : orderedTgroups) 
		{
    			//traversing each class and copying relevant data.
			//key is TranscriptGroup, value is SCTGValue
			//val.barcodeGroup is a nested hash table 
    			bool isKeyPresent = fullEqMap.find_fn(key.first, [&](const SCTGValue& val){
    		      		auto& bg = val.barcodeGroup;
    		      		auto bcIt = bg.find(barcodeIdx);

    		      		// sub-selecting bgroup of this barcode only
    		      		if (bcIt != bg.end()) //if found
				{
    		        		// extracting txp labels
    		        		const std::vector<uint32_t>& txps = key.first.txps;
    		        		txpGroups.emplace_back(txps);
    		        		umiGroups.emplace_back(bcIt->second); //hashmap for this bc
    		        	
					for(auto& ugroup: bcIt->second) //iterate through bc hashmap, adding up counts for each UMI
					{
    		        	  		fragmentCountValidator += ugroup.second;
    		        		}

    		        		// for dumping per-cell eqclass vector
    		        		if(verbose)
					{
    		        			// original counts of the UMI
    		         			uint32_t eqCount {0};
    		          			for(auto& ugroup: bcIt->second)
						{
    		            				eqCount += ugroup.second;
    		        			}
	
	    		         		eqIDs.push_back(static_cast<uint32_t>(key.second)); //transcrip group index/id 
	    		          		eqCounts.push_back(eqCount); //total UMI counts for that Eqc+CB 
	    		        	}
    			      }
		      }); //end lambda 

		      if(!isKeyPresent)
		      {
    		        	jointlog->error("Not able to find key in Cuckoo hash map."
    		    	  	                "Please Report this issue on github");
    		    		jointlog->flush();
    		    		exit(1);
    		      }
    		} //end for

    		if (fragmentCountValidator != umiCounts[barcodeIdx]) 
		{
    			  jointlog->error("Feature count in feature dump doesn't map"
    			                  "with eqclasses frament counts\n");
    			  jointlog->flush();
    			  exit(1);
    		}
		
		
		//Run it
       	 	// perform the UMI deduplication step
       	 	std::vector<SalmonEqClass> salmonEqclasses;
       	 	std::vector<spp::sparse_hash_map<uint16_t, uint16_t>> arboEqClassCount;

       	 	bool dedupOk = dedupClasses(geneAlphas, totalCount, txpGroups,
       	 	                            umiGroups, salmonEqclasses, umiLength,
       	 	                            txpToGeneMap, tiers, gzw, umiEditDistance,
       	 	                            dumpUmiGraph, barcodeStr,
       	 	                            arboEqClassCount,
       	 	                            dumpArborescences,
       	 	                            totalUniEdgesCount, totalBiEdgesCount);
       	 	//std::cout << "\nDedupUMI Run Code: " << dedupOk << "\n";
		if( !dedupOk )
		{
     		   jointlog->error("Deduplication for cell {} failed \n", barcodeStr);
     		   jointlog->flush();
     		   std::exit(74);
     		 }
	
		
		//temp dedup io 
		{
			std::lock_guard<std::mutex> lock(dedupMtx);
			for(auto& it : salmonEqclasses)
			{
				dedupStream << barcodeIdx << "\t";
				dedupStream << it.labels.size() << "\t";

				for(auto& txpId : it.labels)
				{
					dedupStream << txpId << "\t";
				}
				
				dedupStream << it.count << "\n";
			}
			//flush stream here??


		}//end io scope 

		//Progress Update At End
		//printing on screen progress
   		const char RESET_COLOR[] = "\x1b[0m";
   		char green[] = "\x1b[30m";
   		green[3] = '0' + static_cast<char>(fmt::GREEN);
   		char red[] = "\x1b[30m";
   		red[3] = '0' + static_cast<char>(fmt::RED);

   		double cellCount {static_cast<double>(barcode)};
   		if (cellCount > totalCells) 
		{ cellCount = totalCells; }
   		
		double percentCompletion {cellCount*100/totalCells};
   		if (not quiet)
		{
   		  fmt::print(stderr, "\033[A\r\r{}Analyzed {} cells ({}{}%{} of all).{}\n",
        	         green, cellCount, red, round(percentCompletion), green, RESET_COLOR);
   		}
	} //end preprocessing while loop 
}

//jointLog name 


int main(int argc, char* argv[]) 
{
	
	//Driver input params
	#ifdef TAKE_PATHS
	if(argc<2)
        {
                std::cerr << "No arguments provided.  Must provide the path to the file to be read in." << std::endl;
                return 1;
        }
        
	char* bfhPath = argv[1];
        std::cout << "BFH path: " << bfhPath << "\n";
	boost::filesystem::path outputDirectory;

	if(argc>2)
	{
		char* outDir = argv[2];
		//boost::filesystem::path outputDirectory(outDir);
		outputDirectory = outDir;
		std::cout << "Output path: " << outputDirectory.string() << "\n";
	} else {
		outputDirectory = "/mnt/grid/mhammell/hpc/data/data/wunderl/build/methodDev_salmon/test/";
	}
        #endif

	//boost::filesystem::path outputDirectory(argv[2]); //HARDCODE to prevent stupid error 


	//---Begin adding BFH---
	
	/*
	Calling function sig
	int salmonHashQuantify 	( 	AlevinOpts< ProtocolT > &  	aopt,
		bfs::path &  	outputDirectory,
		CFreqMapT &  	freqCounter 
		)
	Original sig
		
		numReads = readBfh( aopt.bfhFile,
                        	    txpNames,
                        	    aopt.protocol.barcodeLength,
                        	    countMap, bcNames,
                        	    freqCounter,
                        	    trueBarcodes,
                        	    totalNormalized,
                        	    hasWhitelist );
	*/
	
	//From Alevin.cpp ~line 900
	CFreqMapT freqCounter;
	freqCounter.reserve(2097152);

	//From AlevinHash.cpp ~line 215
	//NOTE: skipping whitelist code for now!!!
	TrueBcsT trueBarcodes;
	EqMapT fullEqMap;
  	size_t numReads {0};
  	size_t totalNormalized {0};
	std::vector<std::string> txpNames, bcNames;
	uint32_t umiLength = 12; 
	
	{
		std::cout << "\nReading BFH\n";

		//HARDCODED STUFF -- fix later 
		bool hasWhitelist = 0;
		size_t barcodeLength = 16;
	
		std::cout << "UMI length: \t" << umiLength << "\tCB length:\t" << barcodeLength << "\n";

		//set umiLength in static member of Kmer object
		//NOTE: may want to set this inside readBfh instead...
		if (umiLength != alevin::types::AlevinUMIKmer::k()) 
		{
			alevin::types::AlevinUMIKmer::k(umiLength);
      			//aopt.jointLog->info("A custom protocol (END, BC length, UMI length) = ({}, {}, {}) "
                        //   "is being used.  Updating UMI k-mer length accordingly.",
                        //   barEnd, barcodeLength, umiLength);
    		}

		//quick CB test 
		alevin::types::AlevinCellBarcodeKmer::k(barcodeLength);


		boost::filesystem::path bfs_bfhPath {bfhPath};
		

		numReads = readBfh( bfs_bfhPath,
                        	    txpNames,
                        	    barcodeLength,  
                        	    fullEqMap, bcNames,
                        	    freqCounter,
                        	    trueBarcodes,
                        	    totalNormalized,
                        	    hasWhitelist );

		if(numReads == 0)
		{
			std::cerr << "\nError reading BFH\n";
			std::exit(74);
		}

		if(totalNormalized > 0 )
		{
			std::cout << "Skipped " << totalNormalized << " reads due to supplied CB whitelist.\n";
		}

		std::cout << "Detected " << numReads << " total reads in BFH\n";

		#ifdef CHECK_UNIQ_TXP
		{
			spp::sparse_hash_set<std::string> nameHash;
			nameHash.reserve(txpNames.size());
			for(size_t i{0}; i<txpNames.size();i++)
			{
				auto res = nameHash.insert(txpNames[i]);
				if(!res.second)
				{ std::cerr<< "\nWarning: Duplicate annotation found:\t" << txpNames[i] << "\n"; } 

			}
		}
		#endif
		
		#ifdef MY_DEBUG
		printEqMap(fullEqMap);	
		#endif 

	}//end read BFH

	//---End add BFH---
	
	//---Begin preprocessing---
	
	//Alevin reads in t2g 
	/* 
	 * aopt.jointLog->info("Reading transcript to gene Map");
	 *  spp::sparse_hash_map<uint32_t, uint32_t> txpToGeneMap;
	 *  spp::sparse_hash_map<std::string, uint32_t> geneIdxMap;
	 *  getTxpToGeneMap(txpToGeneMap,
	 *                  txpNames,
	 *                  aopt.geneMapFile.string(),
	 *                  geneIdxMap);
	
	 *  GZipWriter gzw(outputDirectory, aopt.jointLog);
	 *  alevinOptimize(bcNames, txpToGeneMap, geneIdxMap,
	 * countMap, aopt, gzw, freqCounter, 0);
	 */

	std::shared_ptr<spdlog::logger> myLogger; //normally aopt.jointLog
        GZipWriter gzw(outputDirectory, myLogger);


	#ifdef TXPMAP_OVERRIDE
	//create dummy txp-gene map; assumes each annotation is a gene 
	spp::sparse_hash_map<uint32_t, uint32_t> txpToGeneMap;
  	spp::sparse_hash_map<std::string, uint32_t> geneIdxMap;	
	

	for (size_t i=0; i<txpNames.size(); i++){
	    txpToGeneMap[i] = i;
	    geneIdxMap[txpNames[i]] = i;
 	}
	#else
	spp::sparse_hash_map<uint32_t, uint32_t> txpToGeneMap;
        spp::sparse_hash_map<std::string, uint32_t> geneIdxMap;

	std::cout << "\nReading in t2g file\n" << std::endl;
	if(argc>3)
	{
		char* t2gPath = argv[3];
		//boost::filesystem::path outputDirectory(outDir);
		std::cout << "T2g path: " << t2gPath << "\n";

		getTxpToGeneMap(txpToGeneMap, txpNames, t2gPath, geneIdxMap);

	} else {
		getTxpToGeneMap(txpToGeneMap, txpNames, "test_t2g.txt", geneIdxMap);
	}

	#endif 
	
	//begin alevinOptimize code from AlevinHash.cpp (there are two methods with the same name!)
	/*Call:
	 * alevinOptimize(bcNames, txpToGeneMap, geneIdxMap,
  			countMap, aopt, gzw, freqCounter, 0);
	 *Sig:
	 *void alevinOptimize( std::vector<std::string>& trueBarcodesVec,
                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                      spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                      EqMapT& fullEqMap,
                      AlevinOpts<ProtocolT>& aopt,
                      GZipWriter& gzw,
                      CFreqMapT& freqCounter,
		      size_t numLowConfidentBarcode)

	*Translation:
	*	bcNames = trueBarcodesVec
	*	countMap = fullEqMap
	*/
	
	std::vector<uint32_t> umiCounts(bcNames.size());
  	for(auto& eq: fullEqMap.lock_table())
	{
     		auto& bg = eq.second.barcodeGroup;
     		for(auto& bcIt: bg)
		{
       			size_t bcCount{0};
       			for(auto& ugIt: bcIt.second)
			{
         			bcCount += ugIt.second;
       			}
       		auto bc = bcIt.first;
       		umiCounts[bc] += bcCount;
  		}
	}

	std::cout << "Starting Preprocessing\n";

	//begin call to optimizer.optimize()
	//This looks like cell filtration based on certain criteria like reads/cell, total riboRNA content, ect
	//Prefiltration before actual processing, we can skip most of this 
	/*
	 * Call:
	 * CollapsedCellOptimizer optimizer;
    	 bool optSuccess = optimizer.optimize(fullEqMap,
    	                                      txpToGeneMap,
    	                                      geneIdxMap,
    	                                      aopt,
    	                                      gzw,
    	                                      trueBarcodesVec,
    	                                      umiCount,
    	                                      freqCounter,
					      numLowConfidentBarcode);
	*
	*Sig:
	*template <typename ProtocolT>
 	 bool CollapsedCellOptimizer::optimize(EqMapT& fullEqMap,
                                       spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                                       spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                                       AlevinOpts<ProtocolT>& aopt,
                                       GZipWriter& gzw,
                                       std::vector<std::string>& trueBarcodes,
                                       std::vector<uint32_t>& umiCount,
				       CFreqMapT& freqCounter,
				       size_t numLowConfidentBarcode)
	*/
	//translation: umiCount=umiCounts
	

	//***Mostly Skipping this for now, not needed for running dedup method****

	size_t totalCells = bcNames.size();
  	size_t totalGenes = geneIdxMap.size();
	std::atomic<uint32_t> curBC{0};
        std::atomic< uint64_t > totalUniEdgesCount(0);
        std::atomic< uint64_t > totalBiEdgesCount(0);

	//CB = 16bp; v2UMI=10bp, v3UMI=12bp
	//Initial Minnow data is v2 
	
	//-HARDCODED for now
  	size_t numWorkerThreads{14}; 
        uint32_t umiEditDistance = 1;
	//moved up to before read bfh: uint32_t umiLength = 10; 
	//-END HARDCODED

	//io options
        bool dumpUmiGraph = 0;
        bool dumpArborescences = 0;
	bool verbose = 1;
	bool quiet = 0;

	
	//Pull the keys (transcrip groups) out of the eq table
	//give each its own unique index value 
 	std::deque<std::pair<TranscriptGroup, uint32_t>> orderedTgroups;
  	//spp::sparse_hash_set<uint64_t> uniqueUmisCounter;
  	uint32_t eqId{0};
  	for(const auto& kv : fullEqMap.lock_table())
	{
  	  // assuming the iteration through lock table is always same
  	  if(kv.first.txps.size() == 1){
  	    orderedTgroups.push_front(std::make_pair(kv.first, eqId));
  	  }
  	  else{
  	    orderedTgroups.push_back(std::make_pair(kv.first, eqId));
  	  }
  	  eqId++;
  	}

	//output header
	boost::filesystem::path dedupOutHeader = outputDirectory / "dedupHeader.txt";	
	std::ofstream headerStream{dedupOutHeader.string()};
	headerStream << txpNames.size() << "\n" << bcNames.size() << "\n";

	for (size_t i=0; i<txpNames.size(); i++)
	{
		headerStream << txpNames[i] << "\n";
        }
	
	for(size_t i=0; i<bcNames.size(); i++)
	{
		headerStream << bcNames[i] << "\n";
	}
	//headerStream.flush();
	headerStream.close();
	
	//---End Preprocessing


	//---Begin Multithreading
	//
	/*void optimizeCell(std::vector<std::string>& trueBarcodes,
                  const std::vector<std::vector<double>>& priorAlphas,
                  std::atomic<uint32_t>& barcode,
                  size_t totalCells, uint32_t umiEditDistance, eqMapT& eqMap,
                  std::deque<std::pair<TranscriptGroup, uint32_t>>& orderedTgroup,
                  std::shared_ptr<spdlog::logger>& jointlog,
                  std::vector<uint32_t>& umiCount,
                  std::vector<CellState>& skippedCB,
                  bool verbose, GZipWriter& gzw, bool noEM, bool useVBEM,
                  bool quiet, std::atomic<double>& totalDedupCounts,
                  std::atomic<uint32_t>& totalExpGeneCounts, double priorWeight,
                  spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                  uint32_t numGenes, uint32_t umiLength, 
                  uint32_t numBootstraps, uint32_t numGibbsSamples,
                  bool naiveEqclass, bool dumpUmiGraph,
                  bool dumpCellEq, bool useAllBootstraps,
                  bool initUniform, CFreqMapT& freqCounter, bool dumpArborescences,
                  spp::sparse_hash_set<uint32_t>& mRnaGenes,
                  spp::sparse_hash_set<uint32_t>& rRnaGenes,
                  std::atomic<uint64_t>& totalUniEdgesCounts,
                  std::atomic<uint64_t>& totalBiEdgesCounts){
	*/
	//
	/*void mtDedupUMIPerCB( EqMapT& fullEqMap,
                      std::vector<std::string>& bcNamess,
                      std::atomic<uint32_t>& barcode,
                      spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                      GZipWriter& gzw,
                      std::shared_ptr<spdlog::logger>& jointlog
                      uint32_t umiLength,
                      uint32_t umiEditDistance,
                      bool dumpUmiGraph,
                      bool dumpArborescences,
                      std::atomic<uint64_t>& totalUniEdgesCounts,
                      std::atomic<uint64_t>& totalBiEdgesCounts)
	*/

	//temp io test
	std::mutex dedupMtx;
	boost::filesystem::path dedupOutPath = outputDirectory / "dedupResults.txt";	
	std::ofstream dedupStream{dedupOutPath.string()};

	std::vector<std::thread> workerThreads;
  	for (size_t tn = 0; tn < numWorkerThreads; ++tn) 
	{
		workerThreads.emplace_back( mtDedupUMIPerCB,
					    std::ref(fullEqMap),
					    std::ref(bcNames),
					    std::ref(curBC),
					    std::ref(orderedTgroups),
					    std::ref(umiCounts),
					    totalCells,
					    totalGenes,
					    std::ref(txpToGeneMap),
					    std::ref(gzw),
					    std::ref(myLogger),
					    umiLength,
					    umiEditDistance,
					    dumpUmiGraph,
					    dumpArborescences,
					    verbose,
					    quiet,
					    std::ref(totalUniEdgesCount),
					    std::ref(totalBiEdgesCount),
					    std::ref(dedupMtx),
					    std::ref(dedupStream));
	}
	
	for (auto& t : workerThreads) 
	{
   		 t.join();
  	}
	
	#if !defined OFF
	//aopt.jointLog->info("Total {0:.2f} UMI after deduplicating", totalDedupCounts);
	std::cout << "Total UMIs after deduplicating:\t" << totalDedupCount <<"\n";
	#endif
	std::cout << "Total bi-directed edges:\t" << totalBiEdgesCount << "\n";
	std::cout << "Total uni-directed edges:\t" << totalUniEdgesCount << "\n";


	std::cout << "\nClosing Streams\n";
	dedupStream.close(); //not 100% necessary 
	gzw.close_all_streams();
	std::cout << "Streams Closed\n";
	//---End Multithreading

	//std::cout << "You found me!!\n";
	
	//DedupClasses caled on line 663 of CollapsedCellOptimizer.cpp
	//Also need initialization of params that starts on line 602
	//Entire loop this is in starts at line 586

	/*
	#if !defined OFF
	//added initiallization
	uint32_t numGenes = 1;
	
	//initialization (start line 603)
	double totalCount{0.0};
     	double totalExpGenes{0};
     	std::vector<uint32_t> eqIDs;
     	std::vector<uint32_t> eqCounts;
     	std::vector<UGroupT> umiGroups;
     	std::vector<tgrouplabelt> txpGroups;
     	std::vector<double> geneAlphas(numGenes, 0.0);
     	std::vector<uint8_t> tiers (numGenes, 0);

	//----Added initialization----
	uint32_t umiLength = 12; //HARDCODED for now, should figure out how to read this option 	
	uint32_t umiEditDistance = 1;
	std::atomic< uint64_t > totalUniEdgesCounts(0);
	std::atomic< uint64_t > totalBiEdgesCounts(0);
	
	//io options
	bool dumpUmiGraph = 0;
	bool dumpArborescences = 0;
	
	std::string trueBarcodeStr = "AAAAGGGGTTTTCCCC"; //need to include string library? -- str lib included as part of alevin utils header 
	
	//https://github.com/greg7mdp/sparsepp#example
	//above--spp::sparse_hash_map<uint32_t, uint32_t > txpToGeneMap;

	//not sure..
	//gzw
	//#include <boost/filesystem.hpp>
	
	//boost::filesystem::path outputDirectory("/mnt/grid/mhammell/hpc/data/data/wunderl/build/methodDev_salmon/test/");
	//boost::filesystem::path outputDirectory(outDir);
	
	//GZipWriter::GZipWriter gzw;
	//above--std::shared_ptr<spdlog::logger> myLogger;
	//aboe--GZipWriter gzw(outputDirectory, myLogger);
		

	//----end-added-----

	//Run it
	// perform the UMI deduplication step
       	std::vector<SalmonEqClass> salmonEqclasses;
       	std::vector<spp::sparse_hash_map<uint16_t, uint16_t>> arboEqClassCount;
       	
	bool dedupOk = dedupClasses(geneAlphas, totalCount, txpGroups,
       	                            umiGroups, salmonEqclasses, umiLength,
       	                            txpToGeneMap, tiers, gzw, umiEditDistance,
       	                            dumpUmiGraph, trueBarcodeStr,
       	                            arboEqClassCount,
       	                            dumpArborescences,
       	                            totalUniEdgesCounts, totalBiEdgesCounts);
	std::cout << "\nDedupUMI Run Code: " << dedupOk << "\n";
	#endif 
	*/
}

