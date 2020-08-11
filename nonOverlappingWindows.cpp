/****************************************************************************
 * nonOverlappingWindows.cpp                                                *
 * Written by Patrick Reilly                                                *
 * Version 1.0 written 2016/07/17                                           *
 * Version 1.1 written 2018/05/15                                           *
 * Version 1.2 written 2018/08/22 Bugfix for header and nonzero omitted stat*
 *   as well as handling a custom statistic column                          *
 * Version 1.3 written 2018/11/08 Filter or use non-N fraction as weight    *
 * Version 1.4 written 2020/01/10 Output sum and denominator if requested   *
 * Description:                                                             *
 *  Calculates the mean of a statistic over non-overlapping windows of      *
 *  user-defined length, and can adjust the denominator of the mean based   *
 *  on a filter column.                                                     *
 *  For instance, when calculating polymorphism in a window, we must        *
 *  exclude sites with Ns from the denominator, and these Ns can be         *
 *  indicated in the filter column by a "1" if present, "0" if absent.      *
 *  Input is a TSV consisting of 3 (4 if a filter column is desired) columns*
 *  identified as scaffold, position, and statistic, respectively.          *
 *  If the scaffold length is not an integral multiple of the window size,  *
 *  the final window's mean is scaled according to its size.                *
 *                                                                          *
 * Syntax: nonOverlappingWindows [options]                                  *
 *  -i:     Path to the input TSV of per-scaffold per-site statistics       *
 *  -o:     Path to desired output TSV of per-scaffold per-window statistics*
 *  -w:     Size of the non-overlapping windows                             *
 *  -n:     Account for a fourth "filter" column to filter certain sites    *
 *          by omitting them from the numerator and denominator.            *
 *  -u:     Output the fraction of usable sites (i.e. not omitted sites) in *
 *          the window as the fourth column.                                *
 *  -s:     Use this column of the input as the statistic to summarize      *
 *          (default: 3, cannot be 1 or 4)                                  *
 *  -f:     Minimum non-N fraction to include in average                    *
 *  -a:     Calculate weighted average using non-N fraction as weight       *
 *  -d:     Output debugging information, including sum and denominator     *
 ****************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <cctype>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <stdexcept>

//Define constants for getopt:
#define no_argument 0
#define required_argument 1
#define optional_argument 2

//Version:
#define version "1.4"

//Usage/help:
#define usage "nonOverlappingWindows\nUsage:\n nonOverlappingWindows [options]\n Options:\n  --input_tsv,-i\tPath to input TSV (default: STDIN)\n  --output_tsv,-o\tPath to output TSV (default: STDOUT)\n  --omit_n,-n\t\tOmit sites indicated in the filter (4th) column\n  --window_size,-w\tSize of the non-overlapping windows\n  --usable_fraction,-u\tOutput the fraction of usable sites in\n\t\t\teach window as column 4\n  --stat_column,-s\tUse this column as the statistic to summarize\n\t\t(default: 3, cannot be 1 or 4)\n  --infimum_nonN,-f\tInfimum fraction of non-Ns to include in average\n\t\t(i.e. include sites with non-N fraction > this value)\n\t\tAssumes column 4 is fraction of non-N bases at site\n  --weighted_average,-a\tCalculate weighted average based on non-N fraction\n\t\tAssumes column 4 is the fraction of non-N bases at the site\n  --debug,-d\t\tOutput debugging information, including sum and denominator as extra columns\n\n Description:\n  Calculates the mean of a statistic over non-overlapping windows\n  across scaffolds in a genome. Sites may be omitted from the average.\n  Input is a 3- or 4-column TSV consisting of scaffold name,\n  position, statistic, and a filter column.\n  If the filter column is 1 and -n is set, the row is omitted from the average.\n  If the fourth column is the fraction of non-N bases,\n  sites may be omitted based on an infimum filter (-f),\n  or a weighted average may be calculated (-a).\n  If the scaffold length is not an integral multiple of the window size,\n  the last window's average is scaled appropriately.\n"

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

string calcDepths(vector<double> &statistics, string scaffold, unsigned long window_size, vector<double> &N_list, bool usable_fraction, unsigned char nonN_weight, double infimum_nonN, bool debug) {
   string output = "";
   if (statistics.size() == 0) { //Fast skip for size 0 scaffolds
      return output;
   }
   if (N_list.size() > 0 && N_list.size() != statistics.size()) { //Something weird happened, we shouldn't ever see this
      throw runtime_error("Somehow there aren't as many weights as statistics for scaffold" + scaffold + ", bugging out.");
   }
   if (debug) {
      cerr << "Processing scaffold " << scaffold << " with omit list of size " << N_list.size() << " and statistic list of size " << statistics.size() << endl;
      if (N_list.size() > 0) {
         if (nonN_weight == 0) {
            cerr << "Omitting sites based on column 4" << endl;
         } else if (nonN_weight == 1) {
            cerr << "Filtering sites unless column 4 >= " << infimum_nonN << endl;
         } else {
            cerr << "Performing weighted average based on column 4" << endl;
         }
      } else {
         cerr << "Nothing in column 4, so performing naive average" << endl;
      }
   }
   double sum = 0.0;
   unsigned long window_start = 1;
   unsigned long last_window;
   double denominator = 0.0;
   unsigned long current_position = 1;
   for (auto stat_iterator = statistics.begin(); stat_iterator != statistics.end(); ++stat_iterator) {
      if (current_position % window_size == 1) {
         window_start = current_position;
         output += scaffold + '\t' + to_string(current_position) + '\t';
      }
      //We either call with empty N_list (to not omit any sites), or N_list.size() == statistics.size()
      if (N_list.size() > 0) {
         if (nonN_weight == 0) {
            denominator += 1.0 - N_list[current_position-1]; //Increment the denominator if omit was 0
            if (N_list[current_position-1] == 0.0) { //If we don't skip this site, add it to the sum
               sum += *stat_iterator;
            }
         } else if (nonN_weight == 1) { //Only include the site if the fraction of non-N bases is high enough, don't include NAs
            if (N_list[current_position-1] > infimum_nonN) {
               sum += *stat_iterator;
               denominator += 1.0;
            }
         } else { //Weight the statistic by the fraction of non-N bases for that site
            sum += *stat_iterator * N_list[current_position-1];
            denominator += N_list[current_position-1];
         }
      } else {
         sum += *stat_iterator;
      }
      if (current_position % window_size == 0) {
         if (N_list.size() == 0) { //No adjustments to denominator if no N_list provided
            denominator = (double)window_size;
         }
         if (denominator > 0.0) { //Avoid dividing by zero
            output += to_string(sum/denominator);
         } else {
            output += "NA";
         }
         if (usable_fraction) {
            output += '\t' + to_string(denominator/window_size);
         }
         if (debug) {
            output += '\t' + to_string(sum) + '\t' + to_string(denominator);
         }
         output += '\n';
         sum = 0.0;
         denominator = 0.0;
      }
      ++current_position;
   }
   //If the scaffold does not contain an integral number of windows:
   if (statistics.size() % window_size > 0) {
      last_window = statistics.size() - window_start + 1;
      if (N_list.size() == 0) { //No adjustments to denominator if no N_list provided
         denominator = (double)last_window;
      }
      if (denominator > 0.0) { //Avoid dividing by zero
         output += to_string(sum/denominator);
      } else {
         output += "NA";
      }
      if (usable_fraction) {
         output += '\t' + to_string(denominator/last_window);
      }
      if (debug) {
         output += '\t' + to_string(sum) + '\t' + to_string(denominator);
      }
      output += '\n';
   }
   return output;
}

int main(int argc, char **argv) {
   //Debug variable:
   bool debug = 0;
   //Variables for processing the TSV:
   bool omit_Ns = 0;
   bool usable_fraction = 0;
   unsigned long int window_size = 10000; //Default window size is 10kb
   unsigned long int stat_column = 3; //Default statistic column is 3
   unsigned char nonN_weight = 0; //Default to treat column 4 as omission indicator
   //1 would be filtering on minimum non-N fraction
   //2 would be weighting average based on non-N fraction
   double infimum_nonN = 0.0; //Default to infimum of 0.0 so we don't include NAs
   string input_tsv = "", output_tsv = "";
   bool use_cin = 1, use_cout = 1; //Default to reading from cin and outputting to cout
   string input_line = "", output_line = "";
   string scaffold_name = "";

   //Variables for getopt_long:
   int optchar;
   int structindex = 0;
   extern int optind;
   //Create the struct used for getopt:
   const struct option longoptions[] {
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'},
      {"debug", no_argument, 0, 'd'},
      {"window_size", required_argument, 0, 'w'},
      {"input_tsv", required_argument, 0, 'i'},
      {"output_tsv", required_argument, 0, 'o'},
      {"omit_n", no_argument, 0, 'n'},
      {"usable_fraction", no_argument, 0, 'u'},
      {"stat_column", required_argument, 0, 's'},
      {"infimum_nonN", required_argument, 0, 'f'},
      {"weighted_average", no_argument, 0, 'a'}
   };
   //Read in the options:
   while ((optchar = getopt_long(argc, argv, "i:o:w:s:nuf:avhd", longoptions, &structindex)) > -1) {
      switch(optchar) {
         case 'v':
            cerr << "nonOverlappingWindows version " << version << endl;
            return 0;
            break;
         case 'h':
            cerr << usage;
            return 0;
            break;
         case 'd':
            cerr << "Debugging mode enabled, outputting sum and denominator columns" << endl;
            debug = 1;
            break;
         case 'i':
            input_tsv = optarg;
            break;
         case 'o':
            output_tsv = optarg;
            break;
         case 'w':
            window_size = atol(optarg);
            break;
         case 'n':
            omit_Ns = 1;
            break;
         case 'u':
            usable_fraction = 1;
            break;
         case 's':
            stat_column = atol(optarg);
            if (stat_column <= 1 || stat_column == 4) {
               cerr << "Statistic column specified as " << optarg << " which is not allowed." << endl;
               cerr << usage;
               return 7;
            }
            break;
         case 'f':
            infimum_nonN = atof(optarg);
            if (infimum_nonN < 0.0 || infimum_nonN >= 1.0) {
               cerr << "Infimum of non-N sites is too small (< 0.0) or too large (>= 1.0), cannot continue." << endl;
               return 8;
            }
            cerr << "Filtering out sites with column 4 less than or equal to " << infimum_nonN << "." << endl;
            nonN_weight = 1;
            break;
         case 'a':
            cerr << "Using column 4 as weight for weighted average." << endl;
            nonN_weight = 2;
            break;
         default:
            cerr << "Unknown option " << (unsigned char)optchar << " supplied." << endl;
            cerr << usage;
            return 1;
            break;
      }
   }

   //Ignore positional arguments
   if (optind < argc) {
      cerr << "Ignoring extra positional arguments starting at " << argv[optind++] << endl;
   }
   
   //Do some error checking on the arguments:
   ifstream input;
   if (input_tsv.length() > 0) { 
      input.open(input_tsv);
      if (!input) {
         cerr << "Error opening input TSV file." << endl;
         return 2;
      }
      use_cin = 0;
   }
   ofstream output;
   if (output_tsv.length() > 0) {
      output.open(output_tsv);
      if (!output) {
         cerr << "Error opening output TSV file." << endl;
         return 3;
      }
      use_cout = 0;
   }
   
   //Set initial state:
   string previous_scaffold = "";
   vector<double> N_list;
   vector<double> scaffold_stats;
   bool header_line = 1;
   unsigned long header_position;
   
   //Now perform the main processing:
   while (getline(use_cin ? cin : input, input_line)) {
      vector<string> line_vector;
      line_vector = splitString(input_line, '\t');
      
      if (header_line) { //Handle a header line if it exists
         header_line = 0; //Only ever check the first line
         try {
            header_position = stoul(line_vector[1]); //Don't do anything if position is numeric
         } catch (const invalid_argument&) { //If position isn't an integer, consider the first line to be a header line
            cerr << "Detected header line and skipped." << endl;
            continue; //Skip the first line
         }
      }
      
      //Test that the TSV is valid and has the necessary columns:
      if (line_vector.size() < 3) { //Too few elements per line
         cerr << "Malformatted input TSV: Less than 3 columns." << endl;
         if (!use_cin) { 
            input.close();
         }
         if (!use_cout) {
            output.close();
         }
         return 4;
      } else if (line_vector.size() < 4 && omit_Ns) { //Need the filter column to perform omission
         cerr << "Malformatted input TSV: Missing filter column." << endl;
         if (!use_cin) {
            input.close();
         }
         if (!use_cout) {
            output.close();
         }
         return 5;
      }
      
      //Convert the column values from strings:
      string scaffold_name = "";
      double local_statistic;
      double omit_position = 0.0;
      scaffold_name = line_vector[0];
      if (stat_column > line_vector.size()) {
         cerr << "Chosen statistic column " << to_string(stat_column) << " is not a valid column in your file." << endl;
         if (!use_cin) {
            input.close();
         }
         if (!use_cout) {
            output.close();
         }
         return 6;
      }
      if (line_vector[stat_column-1] == "NA") {
         local_statistic = 0.0;
         if (nonN_weight > 0) {
            omit_position = 0.0;
         } else {
            omit_position = 1.0;
         }
      } else {
         try {
            local_statistic = stod(line_vector[stat_column-1]);
         } catch (const invalid_argument&) {
            cerr << "Unable to convert stat at " << line_vector[0] << " pos " << line_vector[1] << ": " << line_vector[stat_column-1] << " to double." << endl;
            throw;
         } catch (const out_of_range&) {
            cerr << "Stat at " << line_vector[0] << " pos " << line_vector[1] << ": " << line_vector[stat_column-1] << " is out of range of a double." << endl;
            throw;
         }
         if (line_vector.size() >= 4) {
            try {
               omit_position = stod(line_vector[3]);
            } catch (const invalid_argument&) {
               cerr << "Unable to convert column 4 at " << line_vector[0] << " pos " << line_vector[1] << ": " << line_vector[3] << " to double." << endl;
               throw;
            } catch (const out_of_range&) {
               cerr << "Column 4 at " << line_vector[0] << " pos " << line_vector[1] << ": " << line_vector[3] << " is out of range of a double." << endl;
               throw;
            }
         }
      }
      
      //If we're not on the same scaffold, output the window averages for the previous scaffold:
      if (scaffold_name != previous_scaffold) {
         if (use_cout) {
            cout << calcDepths(scaffold_stats, previous_scaffold, window_size, N_list, usable_fraction, nonN_weight, infimum_nonN, debug);
         } else {
            output << calcDepths(scaffold_stats, previous_scaffold, window_size, N_list, usable_fraction, nonN_weight, infimum_nonN, debug);
         }
         //Make sure to reset the window containers:
         scaffold_stats.clear();
         N_list.clear();
      }
      //Accumulate the window statistics:
      scaffold_stats.push_back(local_statistic);
      if (omit_Ns || nonN_weight) {
         N_list.push_back(omit_position);
      }
      previous_scaffold = scaffold_name;
   }
   //Make sure to capture the last scaffold:
   if (use_cout) {
      cout << calcDepths(scaffold_stats, previous_scaffold, window_size, N_list, usable_fraction, nonN_weight, infimum_nonN, debug);
   } else {
      output << calcDepths(scaffold_stats, previous_scaffold, window_size, N_list, usable_fraction, nonN_weight, infimum_nonN, debug);
   }

   //Close the files if they were used instead of STDIN and STDOUT:
   if (!use_cin) {
      input.close();
   }
   if (!use_cout) {
      output.close();
   }
   
   return 0;
}
