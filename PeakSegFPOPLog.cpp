#include <vector>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
// http://docs.oracle.com/cd/E17076_05/html/programmer_reference/arch_apis.html
#include <dbstl_vector.h>
#include <db_cxx.h>

#include "funPieceListLog.h"

// TODO Use BDB implementation of STL,
// https://blogs.oracle.com/berkeleydb/entry/the_c_standard_template_librar_1
// https://docs.oracle.com/cd/E17276_01/html/programmer_reference/stl.html
// https://docs.oracle.com/cd/E17076_04/html/api_reference/STL/frame_main.html

u_int32_t PiecewiseFunSize(const PiecewisePoissonLossLog&fun){
  return sizeof(PoissonLossPieceLog)*fun.piece_list.size() +
    sizeof(int);
}

void PiecewiseFunCopy(void *dest, const PiecewisePoissonLossLog&fun){
  char *p = (char*)dest;
  int n_pieces = fun.piece_list.size();
  memcpy(p, &n_pieces, sizeof(int));
  p += sizeof(int);
  for(PoissonLossPieceListLog::const_iterator it = fun.piece_list.begin();
      it != fun.piece_list.end(); it++){
    memcpy(p, &(*it), sizeof(PoissonLossPieceLog));
    p += sizeof(PoissonLossPieceLog);
  }
}

void PiecewiseFunRestore(PiecewisePoissonLossLog&fun, const void *src){
  int n_pieces;
  char *p = (char*)src;
  PoissonLossPieceLog piece;
  memcpy(&n_pieces, p, sizeof(int));
  p += sizeof(int);
  for(int piece_i=0; piece_i < n_pieces; piece_i++){
    memcpy(&piece, p, sizeof(PoissonLossPieceLog));
    p += sizeof(PoissonLossPieceLog);
    fun.piece_list.push_back(piece);
  }
}

int main(int argc, char *argv[]){//data_count x 2
  if(argc != 3){
    std::cout << "usage: " << argv[0] << " data.bedGraph penalty\n";
    return 1;
  }
  double penalty = atof(argv[2]);
  std::ifstream bedGraph_file(argv[1]);
  if(!bedGraph_file.is_open()){
    std::cout << "Could not open data file\n";
    return 2;
  }
  std::string line;
  int chromStart, chromEnd, coverage, items, line_i=0;
  double min_log_mean=INFINITY, max_log_mean=-INFINITY, log_data;
  while(std::getline(bedGraph_file, line)){
    line_i++;
    items = sscanf(line.c_str(), "%d %d %d\n", &chromStart, &chromEnd, &coverage);
    if(items!=3){
      printf("error: expected '%%d %%d %%d\\n' on line %d\n", line_i);
      std::cout << line;
      return 3;
    }
    log_data = log(coverage);
    if(log_data < min_log_mean){
      min_log_mean = log_data;
    }
    if(max_log_mean < log_data){
      max_log_mean = log_data;
    }
  }
  int data_count = line_i;
  printf("data_count=%d min_log_mean=%f max_log_mean=%f\n", data_count, min_log_mean, max_log_mean);
  //return 0;
  bedGraph_file.clear();
  bedGraph_file.seekg(0, std::ios::beg);
  dbstl::DbstlElemTraits<PiecewisePoissonLossLog> *funTraits =
    dbstl::DbstlElemTraits<PiecewisePoissonLossLog>::instance();
  funTraits->set_size_function(PiecewiseFunSize);
  funTraits->set_copy_function(PiecewiseFunCopy);
  funTraits->set_restore_function(PiecewiseFunRestore);
  dbstl::db_vector<PiecewisePoissonLossLog> cost_model_mat(data_count * 2);
  //std::vector<PiecewisePoissonLossLog> cost_model_mat(data_count * 2);
  PiecewisePoissonLossLog up_cost, down_cost, up_cost_prev, down_cost_prev;
  PiecewisePoissonLossLog min_prev_cost;
  int verbose=0;
  double cum_weight_i = 0.0, cum_weight_prev_i;
  int data_i = 0;
  double weight;
  while(std::getline(bedGraph_file, line)){
    items = sscanf(line.c_str(), "%d\t%d\t%d\n", &chromStart, &chromEnd, &coverage);
    //if(items != 3)break;
    weight = chromEnd-chromStart;
    //printf("data_i=%d weight=%f coverage=%d\n", data_i, weight, coverage);
    cum_weight_i += weight;
    if(data_i==0){
      // initialization Cdown_1(m)=gamma_1(m)/w_1
      down_cost.piece_list.emplace_back
	(1.0, -coverage, 0.0,
	 min_log_mean, max_log_mean, -1, false);
    }else{
      // if data_i is up, it could have come from down_cost_prev.
      min_prev_cost.set_to_min_less_of(down_cost_prev, verbose);
      int status = min_prev_cost.check_min_of(down_cost_prev, down_cost_prev);
      if(status){
	printf("BAD MIN LESS CHECK data_i=%d status=%d\n", data_i, status);
	printf("=prev down cost\n");
	down_cost_prev.print();
	printf("=min less(prev down cost)\n");
	min_prev_cost.print();
	throw status;
      }
      // C^up_t(m) = (gamma_t + w_{1:t-1} * M^up_t(m))/w_{1:t}, where
      // M^up_t(m) = min{
      //   C^up_{t-1}(m),
      //   C^{<=}_down_{t-1}(m) + lambda/w_{1:t-1}
      // in other words, we need to divide the penalty by the previous cumsum,
      // and add that to the min-less-ified function, before applying the min-env.
      min_prev_cost.set_prev_seg_end(data_i-1);
      min_prev_cost.add(0.0, 0.0, penalty/cum_weight_prev_i);
      // if(data_i==2){
      // 	printf("computing cost data_i=%d\n", data_i);
      // 	verbose=1;
      // }else{
      // 	verbose=0;
      // }
      if(data_i==1){
	up_cost = min_prev_cost;
      }else{
	up_cost.set_to_min_env_of(min_prev_cost, up_cost_prev, verbose);
	status = up_cost.check_min_of(min_prev_cost, up_cost_prev);
	if(status){
	  printf("BAD MIN ENV CHECK data_i=%d status=%d\n", data_i, status);
	  printf("=prev down cost\n");
	  down_cost_prev.print();
	  printf("=min less(prev down cost) + %f\n", penalty);
	  min_prev_cost.print();
	  printf("=prev up cost\n");
	  up_cost_prev.print();
	  printf("=new up cost model\n");
	  up_cost.print();
	  throw status;
	}
      }
      up_cost.multiply(cum_weight_prev_i);
      up_cost.add
	(weight,
	 -coverage*weight,
	 0.0);
      up_cost.multiply(1/cum_weight_i);
      //printf("computing down cost\n");
      // compute down_cost.
      if(data_i==1){
	//for second data point, the cost is only a function of the
	//previous down cost (there is no first up cost).
	down_cost = down_cost_prev;
      }else{
	// if data_i is down, it could have come from up_cost_prev.
	min_prev_cost.set_to_min_more_of(up_cost_prev, verbose);
	status = min_prev_cost.check_min_of(up_cost_prev, up_cost_prev);
	if(status){
	  printf("BAD MIN MORE CHECK data_i=%d status=%d\n", data_i, status);
	  printf("=prev up cost\n");
	  up_cost_prev.print();
	  printf("=min more(prev up cost)\n");
	  min_prev_cost.print();
	  //throw status;
	}
	min_prev_cost.set_prev_seg_end(data_i-1);
	//NO PENALTY FOR DOWN CHANGE
	down_cost.set_to_min_env_of(min_prev_cost, down_cost_prev, verbose);
	status = down_cost.check_min_of(min_prev_cost, down_cost_prev);
	if(status){
	  printf("BAD MIN ENV CHECK data_i=%d status=%d\n", data_i, status);
	  printf("=prev up cost\n");
	  up_cost_prev.print();
	  printf("=min more(prev up cost)\n");
	  min_prev_cost.print();
	  printf("=prev down cost\n");
	  down_cost_prev.print();
	  printf("=new down cost model\n");
	  down_cost.print();
	  throw status;
	}
      }
      down_cost.multiply(cum_weight_prev_i);
      down_cost.add
	(weight,
	 -coverage*weight,
	 0.0);
      down_cost.multiply(1/cum_weight_i);
    }//if(data_i initialization else update
    cum_weight_prev_i = cum_weight_i;
    up_cost_prev = up_cost;
    down_cost_prev = down_cost;
    //printf("data_i=%d data_i+data_count=%d\n", data_i, data_i+data_count);
    cost_model_mat[data_i] = up_cost;
    cost_model_mat[data_i + data_count] = down_cost;
    data_i++;
  }
  //printf("AFTER\n");
  // Decoding the cost_model_vec, and writing to the output matrices.
  double best_cost, best_log_mean, prev_log_mean;
  int prev_seg_end;
  int prev_seg_offset = 0;
  // last segment is down (offset N) so the second to last segment is
  // up (offset 0).
  down_cost = cost_model_mat[data_count*2-1];
  down_cost.Minimize
    (&best_cost, &best_log_mean,
     &prev_seg_end, &prev_log_mean);
  printf("mean=%f end=%d\n", exp(best_log_mean), prev_seg_end);
  // mean_vec[0] = exp(best_log_mean);
  // end_vec[0] = prev_seg_end;
  int out_i=1;
  while(0 <= prev_seg_end){
    // up_cost is actually either an up or down cost.
    up_cost = cost_model_mat[prev_seg_offset + prev_seg_end];
    printf("decoding out_i=%d prev_seg_end=%d prev_seg_offset=%d\n", out_i, prev_seg_end, prev_seg_offset);
    //up_cost.print();
    if(prev_log_mean != INFINITY){
      //equality constraint inactive
      best_log_mean = prev_log_mean;
    }
    up_cost.findMean
      (best_log_mean, &prev_seg_end, &prev_log_mean);
    printf("mean=%f end=%d\n", exp(best_log_mean), prev_seg_end);
    // mean_vec[out_i] = exp(best_log_mean);
    // end_vec[out_i] = prev_seg_end;
    // change prev_seg_offset and out_i for next iteration.
    if(prev_seg_offset==0){
      //up_cost is actually up
      prev_seg_offset = data_count;
    }else{
      //up_cost is actually down
      prev_seg_offset = 0;
    }
    out_i++;
  }//for(data_i
}

