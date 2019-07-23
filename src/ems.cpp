#include "ems.hpp"
#include "ems1.hpp"
#include "ems2.hpp"
#include "ems2p.hpp"
#include <getopt.h>
#include <algorithm>    // std::binary_search, std::sort
#include "utils.h"

void usage(const char *argv) {
	std::cout << "Usage: " << argv << " [OPTIONS] <input-sequence-file>" << std::endl;
	std::cout << "\t-l <l>        Length (l) of (l,d) motif" << std::endl;
	std::cout << "\t-d <d>        Maximum edit distance (d) of (l,d) motif" << std::endl;
	std::cout << "\t-a            Sample rate (n/nPrime)" << std::endl;
	std::cout << "\t-b            AlphabetSize/AlphabetSizePrime" << std::endl;
	exit(-1);
}

int main(int argc, char **argv) {
  std::string input;
  Params params;
  //params.num_threads = omp_get_max_threads();
  params.num_threads = 1;
  params.l = 1;
  params.d = 1;
  double epsilon1 = 1.0, epsilon2 = 1.0;
  int option_char;
  while ((option_char = getopt(argc, argv, "l:d:a:b:")) != -1) {
	switch (option_char)
	  {
		 case 'l': params.l = atoi (optarg); break;
		 case 'd': params.d = atoi (optarg); break;
		 case 'a': epsilon1 = atof (optarg); break;
		 case 'b': epsilon2 = atof (optarg); break;
		 case '?': usage(argv[0]); break;
	  }
  }

  if (argc - optind < 1) {
	usage(argv[0]);
  } 
  input = string(argv[optind]);
  double begin = omp_get_wtime();
  const std::string sample_seqs = removeExtension(input) + "_" + "SampleSeqs" + ".txt";
  const std::string left_seqs = removeExtension(input) + "_" + "LeftSeqs" + ".txt";
  const std::string compress_seqs = removeExtension(input) + "_" + "CompressedSeqs" + ".txt";
  const std::string output_motifs = get_out_file(input, params.l, params.d, "ems2");;
  sample_seq(input, epsilon1, sample_seqs, left_seqs);
  compressInputFile(left_seqs, compress_seqs, epsilon2);
  double end = omp_get_wtime();
  double elapsed1 = end-begin;
  cout << "Sampling and compressing input seqs take " << elapsed1 << " seconds " << endl;
  begin = omp_get_wtime();
  input = sample_seqs;
  std::string output_short = get_out_file(input, params.l, params.d, "ems2");
  Ems2<MotifTreeFast> ems_short(input, params.l, params.d, params);
  ems_short.searchWriteMotifs(params);
  end = omp_get_wtime();
  double elapsed2 = end-begin;
  //cout << "Ems2 on short seqs takes " << elapsed2 << " seconds " << endl;
  begin = omp_get_wtime();
  input = compress_seqs;
  std::string output_compress = get_out_file(input, params.l, params.d, "ems2");
  Ems2<MotifTreeFast> ems_compress(input, params.l, params.d, params);
  ems_compress.searchWriteMotifs(params);
  end = omp_get_wtime();
  double elapsed3 = end-begin;
  //cout << "Ems2 on compressed seqs takes " << elapsed3 << " seconds " << endl;
  begin = omp_get_wtime();
  Reads reads_short, reads_compress, reads_candidate, reads_left, reads_out;
  read_file(output_short.c_str(), reads_short);
  read_file(output_compress.c_str(), reads_compress);
  read_file(left_seqs.c_str(), reads_left);
  std::string sigma = getAlphabet(reads_short);
  for (unsigned int i = 0; i < reads_short.size(); i++) {
    string tmp = reads_short[i];
    compressString(tmp, sigma, epsilon2);
    if (std::binary_search(reads_compress.begin(), reads_compress.end(), tmp)){
        reads_candidate.push_back(reads_short[i]);
        //cout << reads_short[i] << endl;
    }
  }
  end = omp_get_wtime();
  cout << "Candidate motif set size is " << reads_candidate.size() << endl;
  double elapsed4 = end-begin;
  cout << "Merging candidate motifs takes " << elapsed4 << " seconds " << endl;
  begin = omp_get_wtime();
  std::ofstream out(output_motifs);
  for (unsigned int i = 0; i < reads_candidate.size(); i++){
    if (found_in_seqs(reads_candidate[i], reads_left, params.l, params.d)){
        reads_out.push_back(reads_candidate[i]);
        out << reads_candidate[i] << endl;
        //cout << reads_candidate[i] << endl;
    }
  }
  out.close();
  end = omp_get_wtime();
  double elapsed5 = end-begin;
  cout << "Checking candidate motifs takes " << elapsed5 << " seconds " << endl;
  remove(sample_seqs.c_str());
  remove(left_seqs.c_str());
  remove(compress_seqs.c_str());
  remove(output_short.c_str());
  remove(output_compress.c_str());
  double elapsed = elapsed1 + elapsed2 + elapsed3 + elapsed4 + elapsed5;
  cout << "EMS3 found " << reads_out.size() << " motif(s) in "<< elapsed << " seconds " << endl;
}
