#ifndef __PRINT_METHODS__
#define __PRINT_METHODS__

//So far only designed for inclusion as the LAST include in myDriver.cpp


//Experimental
//https://stackoverflow.com/questions/4077609/overloading-output-stream-operator-for-vectort

//Also Consider this: https://stackoverflow.com/a/55270925/1291743

#if 0
template < class T >
std::ostream& operator << (std::ostream& os, const std::vector<T>& v)
{
    os << "[";
    for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    {
        os << " " << *ii <<",";
    }
    os << "]";
    return os;
}
#endif

template < class T >
std::ostream& operator << (std::ostream& os, const std::vector<T>& v)
{
    os << "[";
    for (auto& ele: v)
    {
        os << " " << ele <<",";
    }
    os << "]";
    return os;
}


//Decode UMI hashes
//UMIs are encoded using the pufferfish/Kmer.hpp kmer class
#if 0
std::string toStr() const {
    std::string s(k_, 'X');
    auto& d = data_[0];
    int32_t offset = (2 * k_) - 2;
    for (int32_t idx = 0; offset >= 0; offset -= 2, ++idx) {
      s[idx] = decodeBinary((d >> offset & 0x03));
    }
    return s;
  }
static constexpr char revCodes[4] = {'A', 'C', 'G', 'T'};
#endif

std::string decodeUMI(uint64_t hash)
{
	char revCodes[4] = {'A', 'C', 'G', 'T'};
	
	//hardcoded based on:
	//using AlevinUMIKmer = combinelib::kmers::Kmer<32,2>;
  	//using AlevinCellBarcodeKmer = combinelib::kmers::Kmer<32,5>;
	//note the second param is just a class ID, first param is kmer length (max supported length is 32)
	int32_t k = 32;
	std::string s(k, 'X'); //fill string of length K with Xs
	int32_t offset = (2 * k) - 2;
	for (int32_t idx = 0; offset >= 0; offset -= 2, ++idx) {
      		s[idx] = revCodes[(hash >> offset & 0x03)];
    	}
    	return s;
}

//an attempt at removing the superfluous chars
std::string decodeUMI(uint64_t hash, int32_t k)
{
	char revCodes[4] = {'A', 'C', 'G', 'T'};
	std::string s(k, 'X'); //fill string of length K with Xs
        int32_t offset = (2 * k) - 2;
        for (int32_t idx = 0; offset >= 0; offset -= 2, ++idx) {
                s[idx] = revCodes[(hash >> offset & 0x03)];
        }
        return s;

}




//using  EqMapT = libcuckoo::cuckoohash_map<TranscriptGroup, SCTGValue, TranscriptGroupHasher>
//using  SparseBarcodeMapType = spp::sparse_hash_map< uint32_t, spp::sparse_hash_map< uint64_t, uint32_t > >
// barcodeGroup[barcode][umi]= count;
//using  BarcodeT = uint32_t
//using  UMIT = uint64_t
 
void printEqMap(EqMapT& eqMap)
{
	int32_t umiLen = 10;
        std::cout << "\n";
	for(auto& entry: eqMap.lock_table())
        {
                auto& tg = entry.first;
                auto& sctgval = entry.second;

                //print transcripts
                std::cout << "Transcripts: " << tg.txps <<"\tMass: " << tg.totalMass << "\n";
		std::cout << "-> Count: " << sctgval.count << "\t";
		std::cout << "BarcodeGroup {" << "\n";
		//SparseBarcodeMapType == spp::sparse_hash_map< uint32_t, spp::sparse_hash_map< uint64_t, uint32_t > >
		for(auto& it: sctgval.barcodeGroup) 
		{
			std::cout << "\tBarcode: " << it.first << " { ";
			for(auto& it2: it.second)
			{
				std::cout << "UMI " << decodeUMI(it2.first, umiLen) << " -> Count " << it2.second << "; ";
			}
			std::cout << "}\n";
		}
		std::cout << "}\n\n";

        }
}


#if 0
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
#endif
#endif //End header guard 
