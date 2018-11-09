/****************************************************************************
 * calculateDxy.cpp                                                         *
 * Written by Patrick Reilly                                                *
 * Version 1.0 written 2017/01/29                                           *
 * Version 2.1 written 2017/05/30 (added shared polymorphism and inbred)    *
 * Version 2.2 written 2017/11/13 (no need for list of pseudorefs)          *
 * Version 2.3 written 2018/11/08 (Omit position may output weight instead) *
 *                                                                          *
 * Description:                                                             *
 * This script takes in pseudoreference FASTAs and a TSV describing which   *
 *  pseudoreferences are from which population, and calculates Dxy, Pix,    *
 *  Piy, Dnet, and optionally identifies shared polymorphisms between       *
 *  populations.                                                            *
 * When indicated, all pseudoreferences may be treated as inbred lines, and *
 *  all samples are assumed haploid, where a single allele is chosen at     *
 *  random for each heterozygous site.                                      *
 *                                                                          *
 * Old Syntax: calculateDxy [options] [list of pseudoreference FASTAs]      *
 * Syntax: calculateDxy [options]
 ****************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <cctype>
#include <vector>
#include <map>
#include <sstream>
#include <array>
#include <set>
#include <unordered_map>

//Define constants for getopt:
#define no_argument 0
#define required_argument 1
#define optional_argument 2

//Version:
#define VERSION "2.3"

//Define number of bases:
#define NUM_BASES 4

//Usage/help:
#define USAGE "calculateDxy\nUsage:\n calculateDxy [options]\nOptions:\n -h,--help\tPrint this help\n -v,--version\tPrint the version of this program\n -p,--popfile\tTSV file of FASTA name, and population number\n -s,--shared_poly\tIdentify shared polymorphisms between populations\n -i,--inbred\tTreat pseudoreferences as inbred haploids\n -r,--prng_seed\tSet PRNG seed for random allele selection in inbred lines\n\t\tDefault: 42\n --usable_fraction,-u:\tFourth column represents fraction of unmasked bases\n"

using namespace std;

vector<string> splitString(string line_to_split, char delimiter) {
   vector<string> line_vector;
   string element;
   istringstream line_to_split_stream(line_to_split);
   while (getline(line_to_split_stream, element, delimiter)) {
      line_vector.push_back(element);
   }
   return line_vector;
}

bool openFASTAs(vector<ifstream*> &input_FASTAs, vector<string> input_FASTA_paths) {
   for (auto path_iterator = input_FASTA_paths.begin(); path_iterator != input_FASTA_paths.end(); ++path_iterator) {
      ifstream *input_FASTA = new ifstream(path_iterator->c_str(), ios::in);
      if (!input_FASTA->is_open()) { //Check to make sure the input file was validly opened
         cerr << "Error opening input FASTA: " << *path_iterator << "." << endl;
         return 0;
      }
      input_FASTAs.push_back(input_FASTA);
   }
   return 1;
}

void closeFASTAs(vector<ifstream*> &input_FASTAs) {
   for (auto FASTA_iterator = input_FASTAs.begin(); FASTA_iterator != input_FASTAs.end(); ++FASTA_iterator) {
      (*FASTA_iterator)->close();
      delete *FASTA_iterator;
   }
}

bool readFASTAs(vector<ifstream*> &input_FASTAs, vector<string> &FASTA_lines) {
   unsigned long which_input_FASTA = 0;
   bool ifstream_notfail;
   for (auto FASTA_iterator = input_FASTAs.begin(); FASTA_iterator != input_FASTAs.end(); ++FASTA_iterator) {
      if ((*FASTA_iterator)->fail()) {
         cerr << "Ifstream " << which_input_FASTA+1 << " failed." << endl;
         return 0;
      }
      string FASTA_line;
      ifstream_notfail = (bool)getline(**FASTA_iterator, FASTA_line);
      if (!ifstream_notfail) {
         return ifstream_notfail;
      }
      if (FASTA_lines.size() < input_FASTAs.size()) {
         FASTA_lines.push_back(FASTA_line);
      } else {
         FASTA_lines[which_input_FASTA++] = FASTA_line;
      }
   }
   return ifstream_notfail;
}

string piKey(array<unsigned long, 6> &allele_counts) {
   string pikey = "";
   for (auto array_iterator = allele_counts.begin(); array_iterator != allele_counts.end(); ++array_iterator) {
      pikey += to_string(*array_iterator) + ",";
   }
   return pikey; //We don't really care about the extra comma
}

string dxyKey(array<unsigned long, 6> &pop1_allele_counts, array<unsigned long, 6> &pop2_allele_counts) {
   string dxykey = "";
   for (auto pop1_iterator = pop1_allele_counts.begin(); pop1_iterator != pop1_allele_counts.end(); ++pop1_iterator) {
      dxykey += to_string(*pop1_iterator) + ",";
   }
   for (auto pop2_iterator = pop2_allele_counts.begin(); pop2_iterator != pop2_allele_counts.end(); ++pop2_iterator) {
      dxykey += to_string(*pop2_iterator) + ",";
   }
   return dxykey;
}

bool is_shared_poly(array<unsigned long, 6> &pop1_allele_counts, array<unsigned long, 6> &pop2_allele_counts) {
   unsigned int shared_poly = 0;
   //If any two alleles have non-zero product of frequencies across the two populations, the site is a shared polymorphism
   for (unsigned int i = 0; i < NUM_BASES; i++) {
      shared_poly += ((pop1_allele_counts[i] * pop2_allele_counts[i] > 0) ? 1 : 0);
   }
   return shared_poly > 1;
}

void processScaffold(vector<string> &FASTA_headers, vector<string> &FASTA_sequences, map<unsigned long, unsigned long> &population_map, unsigned long num_populations, unordered_map<string, double> &memoized_pi, unordered_map<string, double> &memoized_dxy, bool shared_poly, bool inbred, bool debug, bool usable) {
   cerr << "Processing scaffold " << FASTA_headers[0].substr(1) << " of length " << FASTA_sequences[0].length() << endl;
   //Do all the processing for this scaffold:
   //Polymorphism estimator: Given base frequencies at site:
   //\hat{\pi} = \(\frac{n}{n-1}\)\sum_{i=1}^{3}\sum_{j=i+1}^{4} 2\hat{p_{i}}\hat{p_{j}}
   //For biallelic sites, this reduces to the standard estimator: \(\frac{n}{n-1}\)2\hat{p}\hat{q}, since
   //p_{3} = p_{4} = 0
   //Dxy estimator is from Nei (1987) Eqn. 10.20 (\hat{d}_{XY} = \Sum_{i,j} \hat{x}_{i} \hat{y}_{j} d_{i,j}
   unsigned long num_sequences = FASTA_sequences.size();
   unsigned long scaffold_length = FASTA_sequences[0].length();
   
   //Containers for various site statistics:
   array<unsigned long, 6> init_base_frequency = { {0, 0, 0, 0, 0, 0} }; //Store the count of A, C, G, T, N, nonN for each site
   vector<array<unsigned long, 6>> population_site_frequencies; //Store population-specific allele counts
   array<double, 4> init_p_hats = { {0.0, 0.0, 0.0, 0.0} }; //Store the estimated allele frequencies for each site
   vector<array<double, 4>> population_p_hats; //Store population-specific estimated allele frequencies
   vector<double> population_pi_hats; //Store population-specific un-corrected polymorphism estimates
   vector<double> D_xys; //Store population pair-specific D_{xy} estimates (absolute, not net divergence)
   vector<double> D_as; //Store population pair-specific D_{a} estimates (net divergence, not absolute)
   vector<bool> SP; //Store population pair-specific indicator of shared polymorphism
   
   for (unsigned long i = 0; i < scaffold_length; i++) {
      //Use site if all populations have at least 2 alleles:
      bool use_site = 1;
   
      //Initialize some of these containers:
      population_site_frequencies.clear();
      population_p_hats.clear();
      population_pi_hats.clear();
      D_xys.clear();
      for (unsigned long i = 0; i < num_populations; i++) {
         population_site_frequencies.push_back(init_base_frequency);
         population_p_hats.push_back(init_p_hats);
         population_pi_hats.push_back(0.0);
      }
      if (debug) {
         cerr << "Counting alleles for site " << i+1 << "." << endl;
      }
      for (unsigned long j = 0; j < num_sequences; j++) {
         switch (FASTA_sequences[j][i]) {
            case 'A':
            case 'a':
               population_site_frequencies[population_map[j]-1][0] += inbred ? 1 : 2; //Add 2 A alleles
               population_site_frequencies[population_map[j]-1][5] += inbred ? 1 : 2; //Add 2 non-N alleles
               break;
            case 'C':
            case 'c':
               population_site_frequencies[population_map[j]-1][1] += inbred ? 1 : 2; //Add 2 C alleles
               population_site_frequencies[population_map[j]-1][5] += inbred ? 1 : 2; //Add 2 non-N alleles
               break;
            case 'G':
            case 'g':
               population_site_frequencies[population_map[j]-1][2] += inbred ? 1 : 2; //Add 2 G alleles
               population_site_frequencies[population_map[j]-1][5] += inbred ? 1 : 2; //Add 2 non-N alleles
               break;
            case 'K': //G/T het site
            case 'k':
               if (inbred) { //Randomly choose one of the alleles
                  population_site_frequencies[population_map[j]-1][rand() <= (RAND_MAX-1)/2 ? 2 : 3]++;
               } else {
                  population_site_frequencies[population_map[j]-1][2]++; //Add 1 G allele
                  population_site_frequencies[population_map[j]-1][3]++; //Add 1 T allele
               }
               population_site_frequencies[population_map[j]-1][5] += inbred ? 1 : 2; //Add 2 non-N alleles
               break;
            case 'M': //A/C het site
            case 'm':
               if (inbred) { //Randomly choose one of the alleles
                  population_site_frequencies[population_map[j]-1][rand() <= (RAND_MAX-1)/2 ? 0 : 1]++;
               } else {
                  population_site_frequencies[population_map[j]-1][0]++; //Add 1 A allele
                  population_site_frequencies[population_map[j]-1][1]++; //Add 1 C allele
               }
               population_site_frequencies[population_map[j]-1][5] += inbred ? 1 : 2; //Add 2 non-N alleles
               break;
            case 'R': //A/G het site
            case 'r':
               if (inbred) { //Randomly choose one of the alleles
                  population_site_frequencies[population_map[j]-1][rand() <= (RAND_MAX-1)/2 ? 0 : 2]++;
               } else {
                  population_site_frequencies[population_map[j]-1][0]++; //Add 1 A allele
                  population_site_frequencies[population_map[j]-1][2]++; //Add 1 G allele
               }
               population_site_frequencies[population_map[j]-1][5] += inbred ? 1 : 2; //Add 2 non-N alleles
               break;
            case 'S': //C/G het site
            case 's':
               if (inbred) { //Randomly choose one of the alleles
                  population_site_frequencies[population_map[j]-1][rand() <= (RAND_MAX-1)/2 ? 1 : 2]++;
               } else {
                  population_site_frequencies[population_map[j]-1][1]++; //Add 1 C allele
                  population_site_frequencies[population_map[j]-1][2]++; //Add 1 G allele
               }
               population_site_frequencies[population_map[j]-1][5] += inbred ? 1 : 2; //Add 2 non-N alleles
               break;
            case 'T':
            case 't':
               population_site_frequencies[population_map[j]-1][3] += inbred ? 1 : 2; //Add 2 T alleles
               population_site_frequencies[population_map[j]-1][5] += inbred ? 1 : 2; //Add 2 non-N alleles
               break;
            case 'W': //A/T het site
            case 'w':
               if (inbred) { //Randomly choose one of the alleles
                  population_site_frequencies[population_map[j]-1][rand() <= (RAND_MAX-1)/2 ? 0 : 3]++;
               } else {
                  population_site_frequencies[population_map[j]-1][0]++; //Add 1 A allele
                  population_site_frequencies[population_map[j]-1][3]++; //Add 1 T allele
               }
               population_site_frequencies[population_map[j]-1][5] += inbred ? 1 : 2; //Add 2 non-N alleles
               break;
            case 'Y': //C/T het site
            case 'y':
               if (inbred) { //Randomly choose one of the alleles
                  population_site_frequencies[population_map[j]-1][rand() <= (RAND_MAX-1)/2 ? 1 : 3]++;
               } else {
                  population_site_frequencies[population_map[j]-1][1]++; //Add 1 C allele
                  population_site_frequencies[population_map[j]-1][3]++; //Add 1 T allele
               }
               population_site_frequencies[population_map[j]-1][5] += inbred ? 1 : 2; //Add 2 non-N alleles
               break;
            case 'N':
               population_site_frequencies[population_map[j]-1][4] += inbred ? 1 : 2; //Add 2 N alleles
               break;
            case '-':
               population_site_frequencies[population_map[j]-1][4] += inbred ? 1 : 2; //Add 2 N alleles
               break;
            default: //Assume that any case not handled here is an N
               population_site_frequencies[population_map[j]-1][4] += inbred ? 1 : 2; //Add 2 N alleles
               break;
         }
      }
      
      //Calculate the total and population-specific allele frequencies:
      if (debug) {
         cerr << "Estimating allele frequencies for site " << i+1 << "." << endl;
      }
      unsigned long population_index = 0;
      unsigned long nonN_bases = 0;
      unsigned long total_bases = 0;
      for (population_index = 0; population_index < num_populations; population_index++) {
         use_site = use_site && population_site_frequencies[population_index][5] >= 2;
         nonN_bases += population_site_frequencies[population_index][5];
         total_bases += population_site_frequencies[population_index][4] + population_site_frequencies[population_index][5];
         string pi_key = piKey(population_site_frequencies[population_index]);
         for (unsigned long j = 0; j < NUM_BASES; j++) {
            population_p_hats[population_index][j] = (double)population_site_frequencies[population_index][j]/(double)population_site_frequencies[population_index][5];
         }
         if (memoized_pi.find(pi_key) != memoized_pi.end()) {
            population_pi_hats[population_index] = memoized_pi[pi_key];
         } else {
            //Calculate \pi_{i} for each population:
            if (debug) {
               cerr << "Estimating pi for each population at site " << i+1 << "." << endl;
            }
            for (unsigned long j = 0; j < NUM_BASES-1; j++) {
               for (unsigned long k = j+1; k < NUM_BASES; k++) {
                  population_pi_hats[population_index] += 2*population_p_hats[population_index][j]*population_p_hats[population_index][k];
               }
            }
            //Do the \frac{n}{n-1} correction, which now makes this Nei (1987) Eqn. 10.5
            if (population_site_frequencies[population_index][5] >= 2) { //Avoid divide-by-zero
               population_pi_hats[population_index] *= (double)population_site_frequencies[population_index][5]/(double)(population_site_frequencies[population_index][5]-1);
            }
            memoized_pi[pi_key] = population_pi_hats[population_index];
         }
      }
      
      //Calculate D_{xy} and D_{a} for each pair of populations, and identify shared polymorphisms:
      if (debug) {
         cerr << "Estimating D_xy and D_a for site " << i+1 << "." << endl;
      }
      for (population_index = 0; population_index < num_populations; population_index++) {
         for (unsigned long population2_index = population_index+1; population2_index < num_populations; population2_index++) {
            string dxy_key = dxyKey(population_site_frequencies[population_index], population_site_frequencies[population2_index]);
            if (memoized_dxy.find(dxy_key) != memoized_dxy.end()) {
               D_xys.push_back(memoized_dxy[dxy_key]);
               double d_net = memoized_dxy[dxy_key] - ((population_pi_hats[population_index] + population_pi_hats[population2_index]) / (double)2.0);
               D_as.push_back(d_net);
            } else {
               double d_xy = 0.0;
               for (unsigned long j = 0; j < NUM_BASES; j++) {
                  for (unsigned long k = 0; k < NUM_BASES; k++) {
                     if (j != k) { //d_{ij} = 1 for i != j, see Nei (1987) Eqn. 10.20
                        d_xy += population_p_hats[population_index][j]*population_p_hats[population2_index][k]; //\hat{x}_{i}\hat{y}_{j} in Nei (1987) Eqn. 10.20
                     }
                  }
               }
               D_xys.push_back(d_xy);
               memoized_dxy[dxy_key] = d_xy;
               double d_net = d_xy - ((population_pi_hats[population_index] + population_pi_hats[population2_index]) / (double)2.0);
               D_as.push_back(d_net);
            }
            //Identify shared polymorphisms:
            if (shared_poly) {
               bool site_is_shared_poly = is_shared_poly(population_site_frequencies[population_index], population_site_frequencies[population2_index]);
               SP.push_back(site_is_shared_poly);
            }
         }
      }

      double usable_fraction = (double)nonN_bases/(double)total_bases;
      
      if (use_site) {
         //Output elements: Scaffold, position, D_{12}, omit site, pi_{i}, D_{ij}, D_{a} values
         cout << FASTA_headers[0].substr(1) << '\t' << i+1;
         population_index = 0;
         //Output D_{12} and omit_site:
         if (!usable) { //If we don't want to output the usable fraction, just output 0
            cout << '\t' << D_xys[0] << '\t' << "0";
         } else {
            cout << '\t' << D_xys[0] << '\t' << usable_fraction;
         }
         //Output \pi_{i} values:
         for (auto population_iterator = population_pi_hats.begin(); population_iterator != population_pi_hats.end(); ++population_iterator) {
            cout << '\t' << *population_iterator;
            population_index++;
         }
         //Output D_{XY}, and D_{a} values:
         unsigned long num_pairs = D_xys.size();
         for (unsigned long pair_index = 0; pair_index < num_pairs; pair_index++) {
            cout << '\t' << D_xys[pair_index];
            cout << '\t' << D_as[pair_index];
         }
         if (shared_poly) {
            for (unsigned long pair_index = 0; pair_index < num_pairs; pair_index++) {
               cout << '\t' << SP[pair_index];
            }
         }
         cout << endl;
      } else { //Do not output any estimators for n < 2
         if (!usable) { //If we don't want to output the usable fraction, just output 1
            cout << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << "0" << '\t' << 1;
         } else {
            cout << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << "0" << '\t' << usable_fraction;
         }
         for (unsigned long j = 1; j <= num_populations; j++) {
            cout << '\t' << "0"; //Output 0 (NA)s for \pi_{i} values as well
         }
         for (unsigned long j = 1; j <= num_populations; j++) {
            for (unsigned long k = j+1; k <= num_populations; k++) {
               cout << '\t' << "0"; //Output 0 (NA)s for D_{XY}
               cout << '\t' << "0"; //Output 0 (NA)s for D_{a}
            }
         }
         if (shared_poly) {
            for (unsigned long j = 1; j <= num_populations; j++) {
               for (unsigned long k = j+1; k <= num_populations; k++) {
                  cout << '\t' << "0"; //Output 0 (NA) for shared polymorphism between i and j
               }
            }
         }
         cout << endl;
      }
   }
}

int main(int argc, char **argv) {
   //Variables for processing the FASTAs:
   vector<string> input_FASTA_paths;
   vector<ifstream*> input_FASTAs;
   
   //Path to file describing which indivs are in which populations:
   string popfile_path;
   
   //Option to identify polymorphisms shared between populations:
   bool shared_poly = 0;
   
   //Handle inbred pseudoreferences as haploids:
   bool inbred = 0;
   unsigned int prng_seed = 42;

   //Option to output debugging info on STDERR
   bool debug = 0;

   //Option to output fraction of usable sites:
   bool usable = 0;
   
   //Variables for getopt_long:
   int optchar;
   int structindex = 0;
   extern int optind;
   //Create the struct used for getopt:
   const struct option longoptions[] {
      {"popfile", required_argument, 0, 'p'},
      {"shared_poly", no_argument, 0, 's'},
      {"inbred", no_argument, 0, 'i'},
      {"prng_seed", required_argument, 0, 'r'},
      {"usable_fraction", no_argument, 0, 'u'},
      {"debug", no_argument, 0, 'd'},
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'}
   };
   //Read in the options:
   while ((optchar = getopt_long(argc, argv, "p:sir:udvh", longoptions, &structindex)) > -1) {
      switch(optchar) {
         case 'p':
            cerr << "Using population TSV file " << optarg << endl;
            popfile_path = optarg;
            break;
         case 's':
            cerr << "Identifying shared polymorphic sites." << endl;
            shared_poly = 1;
            break;
         case 'i':
            cerr << "Assuming all pseudoreferences are haploid." << endl;
            inbred = 1;
            break;
         case 'r':
            cerr << "Setting PRNG seed to " << optarg << endl;
            prng_seed = atoi(optarg);
            break;
         case 'u':
            cerr << "Outputting fraction of usable sites rather than omit column" << endl;
            usable = 1;
            break;
         case 'd':
            cerr << "Outputting debug information." << endl;
            debug = 1;
            break;
         case 'v':
            cerr << "calculateDxy version " << VERSION << endl;
            return 0;
            break;
         case 'h':
            cerr << USAGE;
            return 0;
            break;
         default:
            cerr << "Unknown option " << (unsigned char)optchar << " supplied." << endl;
            cerr << USAGE;
            return 1;
            break;
      }
   }
   
   //Set the seed of the PRNG:
   srand(prng_seed);
   
   //Open the population TSV file:
   ifstream pop_file;
   pop_file.open(popfile_path);
   if (!pop_file) {
      cerr << "Error opening population TSV file " << popfile_path << endl;
      return 9;
   }
   
   //Read the population members into a map from FASTA to population:
   map<unsigned long, unsigned long> population_map;
   set<unsigned long> populations;
   string popline;
   if (debug) {
      cerr << "Reading population TSV file." << endl;
   }
   unsigned long fasta_index = 0;
   while (getline(pop_file, popline)) {
      vector<string> line_vector;
      line_vector = splitString(popline, '\t');
      try {
         population_map[fasta_index++] = stoul(line_vector[1]);
      } catch (const invalid_argument& e) {
         cerr << "Invalid population ID in second column of population TSV." << endl;
         cerr << "Must be a positive integer." << endl;
         pop_file.close();
         return 9;
      }
      input_FASTA_paths.push_back(line_vector[0]);
      populations.insert(stoul(line_vector[1]));
   }
   pop_file.close();
   unsigned long num_populations = populations.size();
   if (debug) {
      cerr << "Read in " << num_populations << " populations." << endl;
   }
   //Output the header line:
   if (!usable) {
      cout << "Scaffold" << '\t' << "Position" << '\t' << "D_1,2" << '\t' << "omit_position";
   } else {
      cout << "Scaffold" << '\t' << "Position" << '\t' << "D_1,2" << '\t' << "site_weight";
   }
   for (unsigned long i = 1; i <= num_populations; i++) {
      cout << '\t' << "pi_" << i;
   }
   for (unsigned long i = 1; i <= num_populations; i++) {
      for (unsigned long j = i+1; j <= num_populations; j++) {
         cout << '\t' << "D_" << i << ',' << j;
         cout << '\t' << "Da_" << i << ',' << j;
      }
   }
   if (shared_poly) {
      for (unsigned long i = 1; i <= num_populations; i++) {
         for (unsigned long j = i+1; j <= num_populations; j++) {
            cout << '\t' << "Shared_Poly_" << i << ',' << j;
         }
      }
   }
   cout << endl;
   
   //Open the input FASTAs:
   bool successfully_opened = openFASTAs(input_FASTAs, input_FASTA_paths);
   if (!successfully_opened) {
      closeFASTAs(input_FASTAs);
      cerr << "Unable to open at least one of the FASTAs provided." << endl;
      return 2;
   }
   cerr << "Opened " << input_FASTAs.size() << " input FASTA files out of " << input_FASTA_paths.size() << " paths provided." << endl;
   
   //Set up the vector to contain each line from the n FASTA files:
   vector<string> FASTA_lines;
   FASTA_lines.reserve(input_FASTA_paths.size());
   
   //Set up maps to memoize pi and Dxy:
   unordered_map<string, double> memoized_pi;
   unordered_map<string, double> memoized_dxy;

   //Iterate over all of the FASTAs synchronously:
   vector<string> FASTA_headers;
   FASTA_headers.reserve(input_FASTA_paths.size());
   vector<string> FASTA_sequences;
   FASTA_sequences.reserve(input_FASTA_paths.size());
   while (readFASTAs(input_FASTAs, FASTA_lines)) {
      //Check if we're on a header line:
      bool all_header_lines = 1;
      bool any_header_lines = 0;
      for (auto line_iterator = FASTA_lines.begin(); line_iterator != FASTA_lines.end(); ++line_iterator) {
         all_header_lines = all_header_lines && ((*line_iterator)[0] == '>');
         any_header_lines = any_header_lines || ((*line_iterator)[0] == '>');
      }
      if (all_header_lines) {
         if (!FASTA_sequences.empty()) {
            processScaffold(FASTA_headers, FASTA_sequences, population_map, num_populations, memoized_pi, memoized_dxy, shared_poly, inbred, debug, usable);
            FASTA_sequences.clear();
         }
         FASTA_headers = FASTA_lines;
         for (auto header_iterator = FASTA_headers.begin()+1; header_iterator != FASTA_headers.end(); ++header_iterator) {
            if (*header_iterator != FASTA_headers[0]) {
               cerr << "Error: FASTAs are not synchronized, headers differ." << endl;
               cerr << FASTA_headers[0] << endl;
               cerr << *header_iterator << endl;
               closeFASTAs(input_FASTAs);
               return 3;
            }
         }
      } else if (any_header_lines) {
         cerr << "Error: FASTAs are not synchronized, or not wrapped at the same length." << endl;
         closeFASTAs(input_FASTAs);
         return 4;
      } else {
         unsigned long sequence_index = 0;
         for (auto line_iterator = FASTA_lines.begin(); line_iterator != FASTA_lines.end(); ++line_iterator) {
            if (FASTA_sequences.size() < FASTA_lines.size()) {
               FASTA_sequences.push_back(*line_iterator);
            } else {
               FASTA_sequences[sequence_index++] += *line_iterator;
            }
         }
      }
   }
   
   //Catch any IO errors that kicked us out of the while loop:
   unsigned long which_input_FASTA = 0;
   for (auto FASTA_iterator = input_FASTAs.begin(); FASTA_iterator != input_FASTAs.end(); ++FASTA_iterator) {
      if ((*FASTA_iterator)->bad()) {
         cerr << "Error reading input FASTA: " << input_FASTA_paths[which_input_FASTA] << endl;
         cerr << "Fail bit: " << ((*FASTA_iterator)->rdstate() & ifstream::failbit) << " Bad bit: " << (*FASTA_iterator)->bad() << " EOF bit: " << (*FASTA_iterator)->eof() << endl;
         closeFASTAs(input_FASTAs);
         return 5;
      }
      which_input_FASTA++;
   }
   //If no errors kicked us out of the while loop, process the last scaffold:
   processScaffold(FASTA_headers, FASTA_sequences, population_map, num_populations, memoized_pi, memoized_dxy, shared_poly, inbred, debug, usable);
   
   //Close the input FASTAs:
   closeFASTAs(input_FASTAs);
   
   return 0;
}
