#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// filter triple that has monotonic pattern and valley pattern
NumericVector filt_pattn(NumericVector peaks, NumericVector stat, int dist) {
  int len = peaks.size();
  int i=1, a, b, c;

  //remove valley pattern
  while(i < (len-1)){
    a = peaks[i] - 1;
    b = peaks[i - 1] - 1;
    c = peaks[i + 1] - 1;
    if((c - b <= dist) && (stat[a] < stat[b]) && (stat[a] < stat[c])){
      IntegerVector idx = seq_len(len) - 1;
      // Rcout << i << std::endl;
      peaks = peaks[idx != i];
      len = peaks.size();
      i = 1;
    }
    else{
      i += 1;
    }
  }

  // remove monotonic pattern
  i = 1;
  while(i < (len - 1)){
    a = peaks[i] - 1;
    b = peaks[i - 1] - 1;
    c = peaks[i + 1] - 1;
    if((a - b <= dist) && (c - a <= dist) && (stat[a] > stat[b]) && (stat[c] > stat[a])){
        IntegerVector idx = seq_len(len) - 1;
        // Rcout << i << std::endl;
        peaks = peaks[idx != i];
        len = peaks.size();
        i = 1;
    }
    else{
      i += 1;
    }
  }
  return peaks;
}


// [[Rcpp::export]]
NumericVector combPeaks(NumericVector peaks, NumericVector stat, int dist) {
  // recursively combine peaks that are too close

  // filter out pattern
  peaks = filt_pattn(peaks, stat, dist);

  int nn = peaks.size();
  int current_peak, next_peak, i=0;
  std::vector<int> new_peaks;
  std::vector<int> del_peaks;

  current_peak = peaks[0];
  new_peaks.push_back(current_peak);

  while(i < nn-1){
    next_peak = peaks[i + 1];
    // whether current peak and next peak are close
    if(next_peak - current_peak <= dist){
     // which one has bigger stat
     if(stat[current_peak] > stat[next_peak]){
       del_peaks.push_back(next_peak);
       i++;
     }else{
       del_peaks.push_back(current_peak);
       new_peaks.pop_back();
       new_peaks.push_back(next_peak);
       current_peak = next_peak;
       i++;
     }
    }else{
      new_peaks.push_back(next_peak);
      current_peak = next_peak;
      i++;
    }
  }

  return wrap(new_peaks);
}


