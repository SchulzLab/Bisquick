// ==========================================================================
//                                  bisquick
// ==========================================================================
//
// ==========================================================================
// Author: Luis Enrique Ramirez Chavez <lramirez@mpi-inf.mpg.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/seq_io.h>
#include <seqan/alignment_free.h>
#include <seqan/arg_parse.h>

#include <dirent.h>
#include <iostream>
#include <limits>
#include <thread>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <utility>
#include <stack>
#include <algorithm>

using namespace seqan;

static const char gendir[] = "/MMCI/MS/MethylationAlgos/work/BisulfiteAnalysis/lab/genome3";
//static const char gendir[] = "/MMCI/MS/MethylationAlgos/work/BisulfiteAnalysis/dev/seqan-trunk-build/debug/bin/genoma/";

// ==========================================================================
// GLOBAL VARIABLES
// ==========================================================================

// ksize is the size of the kmer
unsigned int ksize = 25;
// number of cgs
unsigned int ncpgs;
// apparently dummy 
unsigned char *kmers;



// genome vector that contain the genome in each element a chromosome
std::vector< String<Dna> >genome;

// Intervals per  chromosome
std::vector<uint>interpchrom;


// cg_arrays: contains 4 vectors to map the cg in the genome
// chromosome_nr: contains the number of chromosome that 
struct cg_arrays {
  std::vector<unsigned char>chromosome_nr;
  std::vector<unsigned int>cg_position_inchr;
  std::vector<unsigned int>methylated;
  std::vector<unsigned int>unmethylated;
} cg_index;


// intervals_index is the index to map all the intervals
struct intervals_arrays {
  std::vector<unsigned int>start_interval;
  std::vector<unsigned int>end_interval;
  std::vector<unsigned int>cg_start;
  std::vector<unsigned int>cg_offset;
  //  std::vector< std::pair <unsigned int,unsigned int> > intervals;
  //std::vector< std::pair <unsigned int,unsigned int> > cg_ref;
} intervals_index;

//map intervals and objects
struct map_value {
  std::vector<unsigned int>interval_index;
  std::vector<unsigned int>offset;
} mapval;

//comprss map
std::map<Dna5String, map_value> commap;



// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct AppOptions {
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;
    unsigned ksize;
    // The first (and only) argument of the program is stored here.
    seqan::CharString genomedir;
    seqan::CharString readsdir;
    seqan::CharString output;
    AppOptions() :
     verbosity(1),
     ksize(5),
     genomedir(gendir)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

/*FUNCTION DECLARATIONS*/

void fill_indices(seqan::CharString seq, int fn);
int readfasta(char * fastafile, int fn);
void read_files(char * curdir);
void read_files(char * curdir,int dirtype);
void open_directory(char * curdir);

Dna5String kmer_compress(Dna5String original);
int  openreads();
void create_maps();
void print_mainstructs();
void prubea(int a);
void count_methylation(Dna5String ref,  Dna5String currentkmer,int cgnum ,  int intum, int kmerpos);
void print_mainstructs();



// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv) {
  // Setup ArgumentParser.
  seqan::ArgumentParser parser("bisquick");
  // Set short description, version, and date.
  setShortDescription(parser, "Bisulfite analysis as fast as baking a pie!");
  setVersion(parser, "0.1");
  setDate(parser, "December 2014");
  // Define usage line and long description.
  addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
  addDescription(parser, 
    "Bisquick: Efficient algorithm methylation rate calculation.");
  // We require one argument.
  addOption(parser, 
    seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
  addOption(parser, 
    seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
  addOption(parser, 
    seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));
  addOption(parser, 
    seqan::ArgParseOption(
      "k", "ksize", "Length of the kmer.", 
      seqan::ArgParseArgument::INTEGER, "INT"));
  addOption(parser, seqan::ArgParseOption(
    "g", "genomedir", "Directory of the genome.",
    seqan::ArgParseArgument::STRING, "STRING"));
  addOption(parser, seqan::ArgParseOption(
    "r", "readsdir", "Directory of the reads.",
    seqan::ArgParseArgument::STRING, "STRING"));
  addOption(parser, seqan::ArgParseOption(
    "o", "output", "Output file.",
    seqan::ArgParseArgument::STRING, "STRING"));
  setDefaultValue(parser, "ksize", "25");

  // Add Examples Section.
  addTextSection(parser, "Examples");
  addListItem(parser, "\\fBbisquick\\fP \\fB-v\\fP \\fItext\\fP",
              "Call with \\fITEXT\\fP set to \"text\" with verbose output.");
  // Parse command line.
  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

  // Only extract  options if the program will continue after parseCommandLine()
  if (res != seqan::ArgumentParser::PARSE_OK)
      return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValue(options.ksize, parser, "ksize");
    getOptionValue(options.genomedir, parser, "genomedir");
    getOptionValue(options.readsdir, parser, "readsdir");
    getOptionValue(options.output, parser, "output");
    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function open_directory(char * curdir);
// --------------------------------------------------------------------------

// receives a  directory curdir to open all files in it */
void open_directory(char * curdir)
{

  DIR           *d;
  struct dirent *dir;
  d = opendir(curdir);
  int file_counter=0;			 // file counter
  std::cout << "___________________"<<curdir<<"\n";
  if (d){
  std::cout << "__arf____________________________________________________________________\n";
    while ((dir = readdir(d)) != NULL){
      if(dir->d_name)
	if((std::strcmp(dir->d_name,".")) && (std::strcmp(dir->d_name,"..")) ) {
	  int totsize = strlen(curdir)+strlen(dir->d_name);
	  char * fastafile=NULL;
	  fastafile = (char *)malloc(sizeof(char) * totsize);
	  strcpy (fastafile,curdir);
	  strcat (fastafile,"/");
	  strcat (fastafile,dir->d_name);
	  //if(filetype == 1)
	  //openreads();
	  readfasta(fastafile,file_counter);
	  free(fastafile);
	  file_counter++;
	}
    }
    closedir(d);
  }
  else{
    std::cerr<<"Invalid directory"<<std::endl;
  }

}


// --------------------------------------------------------------------------
// Function  create_maps()
// --------------------------------------------------------------------------
// function to create maps index


void create_maps(){

  for(uint j = 0; j != intervals_index.start_interval.size(); j++) {
    uint offs = 0;
    std::cout<<"int: "<<j<<" of "<<intervals_index.start_interval.size()<<"\n";
    std::cout<<"intervalos: "<<intervals_index.end_interval[j]<<" ksize "<<ksize<<std::endl;
    for(uint i = intervals_index.start_interval[j]; 
	i<=intervals_index.end_interval[j]-ksize; i++) {
	Dna5String currentkmer;
	Dna5String currcompresed;
	uint curr_chrm =(uint)cg_index.chromosome_nr[intervals_index.cg_start[j]]; 
	currentkmer = infix(genome[curr_chrm],i,i+ksize);
	currcompresed = currentkmer;
	int Tcounter = 0;
	//std::cout<<"aqui: "<<i<<" ~ "<<intervals_index.end_interval[j]-ksize<<std::endl;
	for (unsigned k = 0; k < ksize; ++k) {
	  if(currcompresed[k] == (Dna5)'T') Tcounter++;
	  if(currcompresed[k] == (Dna5)'C'){
	    currcompresed[k] = (Dna5)'T';
	    Tcounter++;	
	  }
	}
	typedef std::map<Dna5String, map_value> MapType;    // Your map type may vary, just change the typedef
	Dna5String k = currcompresed;   // assume we're searching for keys equal to 4
	map_value v  = {{j},{offs}};
	MapType::iterator lb = commap.lower_bound(k);
	if(lb != commap.cend() && !(commap.key_comp()(k, lb->first))) {
	  lb->second.interval_index.emplace_back(j);
	  lb->second.offset.emplace_back(offs);
	  // key already exists
	  // update lb->second if you care to
	}
	else {
	  // the key does not exist in the map
	  // add it to the map
	  commap.insert(lb,  MapType::value_type(k, v));
	  //std::pair<Dna5String, map_value>(currcompresed,varm));
	  // mymap.insert(lb, MapType::value_type(k, v));    // Use lb as a hint to insert,
                                                    // so it can avoid another lookup
	}
	//commap.insert(currcompresed, {{2},{2}});
		      //	mymap.insert(pair<int, vector<int> > (10, {1, 2, 3}));
		      //	commap[currcompresed]={}
	offs++;
    }
  }
}




// --------------------------------------------------------------------------
// Function  print_mainstructs()
// --------------------------------------------------------------------------
// for debug 
#define VAR_NAME(n) #n
template <class T> 
void printvec(std::vector<T> const &myvec, const std::string &vecname ) {
  typename std::vector<T>::const_iterator it;
  std::cout <<vecname<< " vector contains:";
  if (typeid(unsigned char) == typeid(T)) 
    for( it = myvec.begin(); it < myvec.end(); it++) {
      std::cout << " " <<(int) *it;
    } else {
    for( it = myvec.begin(); it < myvec.end(); it++) {
      std::cout << " " << *it;
    }
  }
  std::cout << std::endl;
}


void print_mainstructs(){
  std::cout << std::string(30, '*') << std::endl;
  std::cout<<cg_index.chromosome_nr.size()<<std::endl;
  printvec(cg_index.chromosome_nr,VAR_NAME(cg_index.chromosome_nr));
  printvec(cg_index.cg_position_inchr,VAR_NAME(cg_index.cg_position_inchr));
  printvec(cg_index.methylated,VAR_NAME(cg_index.methylated));
  printvec(cg_index.unmethylated,VAR_NAME(cg_index.unmethylated));
  printvec(intervals_index.start_interval,VAR_NAME(intervals_index.start_interval));
  printvec(intervals_index.end_interval,VAR_NAME(intervals_index.end_interval));
  printvec(intervals_index.cg_start,VAR_NAME(intervals_index.cg_start));
  printvec(intervals_index.cg_offset,VAR_NAME(intervals_index.cg_offset));
  std::cout<<"genome size: "<<genome.size()<<std::endl;
  std::cout<<"genome 1: "<<genome[0]<<std::endl;
  std::cout<<"genome 2: "<<genome[1]<<std::endl;
  std::cout<<"genome 1: "<<genome[2]<<std::endl;
  std::cout<<"genome 2: "<<genome[3]<<std::endl;

  std::cout << std::string(30, '*') << std::endl;
  // for(auto it = commap.cbegin(); it != commap.cend(); ++it) {
  //   std::cout << it->first << "---> {";
  //   for(uint j = 0; j<it->second.interval_index.size(); j++)
  //     std::cout<<"{" <<it->second.interval_index[j] << "," << it->second.offset[j] << "}";
  //   std::cout<<"}\n";
  // }
}

// --------------------------------------------------------------------------
// Function fill_indices(seq, fn)
// --------------------------------------------------------------------------
// this function fills cg_index, interval_index and interpchrom



void fill_indices(seqan::CharString seq, int fn){
  int newfile = 1;
  seqan::CharString id;



  genome.emplace_back(seq);
  std::string s;
  std::cout << id <<std::endl;
  CharString needle = "CG";
  Finder<CharString> finder(seq);
  Pattern<CharString,  ShiftAnd> pattern(needle);
  int newinterval = 0;
  // std ::cout<<"tamano de seq:"<< length(seq)<<std ::endl;
  while (find(finder, pattern)) {
    ncpgs++;
    uint cpos = beginPosition(finder);
    uint gpos = endPosition(finder);
    // save the position and the chromosome where "CG" is found
    cg_index.cg_position_inchr.emplace_back(cpos);
    cg_index.chromosome_nr.emplace_back(fn);
    //L directory.emplace_back(std::make_pair(fn,beginPosition(finder)));
    int ini = (int)(cpos)-ksize+2;
    uint fin = gpos+ksize-2;

    if (ini < 0) {
	ini = 0;
	fin = ksize-1;
    }
    if (fin > length(seq)) {
      fin = ksize-1;
      ini =length(seq)-ksize;
    }

    if(intervals_index.start_interval.size()>0) {
      if(newfile) {
	intervals_index.start_interval.emplace_back(ini);
	intervals_index.end_interval.emplace_back(fin);
	intervals_index.cg_start.emplace_back(cg_index
					      .cg_position_inchr.size()-1);
	intervals_index.cg_offset
	               .emplace_back(1);
	newinterval = 0;
	newfile = 0; 
      } else { 
	  if(intervals_index.end_interval.back() >= gpos) {
	    intervals_index.end_interval.back() = gpos+ksize-2;
	    newinterval++;
	    intervals_index.cg_offset.back()=newinterval+1;
	  } else {
	    intervals_index.start_interval.emplace_back(ini);
	    intervals_index.end_interval.emplace_back(fin);
	    intervals_index.cg_start.emplace_back(cg_index.cg_position_inchr.size()-1);
	    intervals_index.cg_offset.emplace_back(1);
	    newinterval = 0;
	  }
      }
    } else {
      intervals_index.start_interval.emplace_back(ini);
      intervals_index.end_interval.emplace_back(fin);
      intervals_index.cg_start.emplace_back(cg_index.cg_position_inchr.size()-1);
      intervals_index.cg_offset.emplace_back(1);
      newinterval = 0;
      newfile = 0;
    }
  }
  cg_index.methylated.resize(cg_index.cg_position_inchr.size());
  cg_index.unmethylated.resize(cg_index.cg_position_inchr.size());
  interpchrom.emplace_back(intervals_index.start_interval.size()-1);

}


// --------------------------------------------------------------------------
// Function readfasta(fastafile, fn)
// --------------------------------------------------------------------------
// reads a file fastafile in fasta format and the file number 
// fn to identify it */

int readfasta(char * fastafile, int fn) {
  kmers = nullptr;
  int newfile = 1;
  seqan::CharString id;
  //String<Dna> seq;
  std::cout <<"fastaf2: " <<fastafile << ::std::endl;
  seqan::CharString seq;
  seqan::SequenceStream seqStream(fastafile);
  if (!isGood(seqStream)) {
    std::cerr << "ERROR: Could not open the file "<<fastafile<<"!.\n";
    return 1;
  }
  if (readRecord(id, seq, seqStream) != 0) {
    std::cerr << "ERROR: Could not read from "<<fastafile<<"!\n";
    return 1;
  }
  // Save the chromosome
  fill_indices(seq,fn*2);
  String<Dna> revseq = seq;
  seqan::CharString seq2rev = DnaStringReverseComplement(revseq);
  fill_indices(seq2rev,fn*2+1);

  // std::cout << myString << ::std::endl;
  // std::cout << DnaStringReverseComplement(myString) << ::std::endl;
  // seqan::CharString seq2 = myString;
  // seqan::CharString seq2rev = DnaStringReverseComplement(myString);
  // String<Dna> myStringr = seq2rev;

  // std::cout << "Resultados:" << ::std::endl;
  // std::cout << myString << ::std::endl;
  // std::cout << myStringr << ::std::endl;
  // std::cout << seq2 << ::std::endl;
  // std::cout << seq2rev << ::std::endl;



  return 0;
  genome.emplace_back(seq);


  std::string s;
  std::cout << id <<std::endl;
  CharString needle = "CG";
  Finder<CharString> finder(seq);
  Pattern<CharString,  ShiftAnd> pattern(needle);
  int newinterval = 0;
  while (find(finder, pattern)) {
    ncpgs++;
    uint cpos = beginPosition(finder);
    uint gpos = endPosition(finder);
    // save the position and the chromosome where "CG" is found
    cg_index.cg_position_inchr.emplace_back(cpos);
    cg_index.chromosome_nr.emplace_back(fn);
    //L directory.emplace_back(std::make_pair(fn,beginPosition(finder)));
    if(intervals_index.start_interval.size()>0) {
      if(newfile) {
	intervals_index.start_interval.emplace_back((int)(cpos)-ksize+2);
	intervals_index.end_interval.emplace_back(gpos+ksize-2);
	intervals_index.cg_start.emplace_back(cg_index.cg_position_inchr.size()-1);
	intervals_index.cg_offset.emplace_back(1);
	newinterval = 0;
	newfile = 0; 
      } else { 
	  if(intervals_index.end_interval.back() >= gpos) {
	    intervals_index.end_interval.back() = gpos+ksize-2;
	    newinterval++;
	    intervals_index.cg_offset.back()=newinterval+1;
	  } else {
	    intervals_index.start_interval.emplace_back(cpos-ksize+2);
	    intervals_index.end_interval.emplace_back(gpos+ksize-2);
	    intervals_index.cg_start.emplace_back(cg_index.cg_position_inchr.size()-1);
	    intervals_index.cg_offset.emplace_back(1);
	    newinterval = 0;
	  }
      }
    } else {
      int ini = (int)(cpos)-ksize+2;
      if (ini < 0) 
	ini = 0;
      intervals_index.start_interval.emplace_back(ini);
      intervals_index.end_interval.emplace_back(gpos+ksize-2);
      intervals_index.cg_start.emplace_back(cg_index.cg_position_inchr.size()-1);
      intervals_index.cg_offset.emplace_back(1);
      newinterval = 0;
      newfile = 0;
    }
  }
  cg_index.methylated.resize(cg_index.cg_position_inchr.size());
  cg_index.unmethylated.resize(cg_index.cg_position_inchr.size());
  interpchrom.emplace_back(intervals_index.start_interval.size()-1);
  return 0;
}





// --------------------------------------------------------------------------
// Function main()
/// --------------------------------------------------------------------------
// Program entry point.

int main(int argc, char const ** argv) {
  // Parse the command line.
  seqan::ArgumentParser parser;
  AppOptions options;
  seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
  // If there was an error parsing or built-in argument parser functionality
  // was triggered then we exit the program.  The return code is 1 if there
  // were errors and 0 if there were none.
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;
  // Print the command line arguments back to the user.
  if (options.verbosity > 0) {
    std::cout << "__OPTIONS____________________________________________________________________\n"
              << '\n'
	      << "output file     \t" << options.output << '\n'
              << "VERBOSITY\t" << options.verbosity << '\n'
              << "ksize   \t" << options.ksize << '\n'
              << "Directory of the genome: \t" << options.genomedir << '\n'
              << "Directory of the reads: \t" << options.readsdir << "\n"<<std::flush;
  }
    

  ksize = options.ksize;
  ncpgs = 0;

  open_directory(toCString(options.genomedir));
  
  //  print_mainstructs();

  create_maps();

  // std::cout<<"aaaaaagenome size: "<<genome.size()<<std::endl;
  // //  std::cout<<"genome 1: "<<toCString(genome[0])<<std::endl;
  // std::cout<<"genome 1: "<<genome[1]<<std::endl;
  // //  std::cout<<"genome 2: "<<reverseComplement(genome[1])<<std::endl;
  // std::cout << std::string(30, '*') << std::endl; 


  // String<Dna> myString = "attacgg";
  // std::cout << myString << ::std::endl;
  // std::cout << DnaStringReverseComplement(myString) << ::std::endl;
  // seqan::CharString seq2 = myString;
  // seqan::CharString seq2rev = DnaStringReverseComplement(myString);
  // String<Dna> myStringr = seq2rev;

  // std::cout << "Resultados:" << ::std::endl;
  // std::cout << myString << ::std::endl;
  // std::cout << myStringr << ::std::endl;
  // std::cout << seq2 << ::std::endl;
  // std::cout << seq2rev << ::std::endl;
 
  // print_mainstructs();
    return 0;
}
