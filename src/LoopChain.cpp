//Contains functions related to the use of loop chaining for interloop data optimization

#include "RAJA/pattern/forall.hpp"
#include "RAJA/util/LoopChain.hpp"
#include <vector>
#include <string>
#include <map>
namespace RAJA
    {


  isl_stat print_point_val(isl_point *pnt, void *user) {
    isl_val * pointVal = isl_point_get_coordinate_val(pnt, isl_dim_set, 0);
    isl_printer * p = isl_printer_to_file((isl_ctx*) user, stdout);
    std::cout << "pointVal: ";
    isl_printer_print_val(p,pointVal);
    return isl_stat_ok;
  }
  isl_stat collect_point_val(isl_point *pnt, void *user) {
    
    std::vector<int> * container = (std::vector<int>*) user;
    isl_val * pointVal = isl_point_get_coordinate_val(pnt, isl_dim_set, 0);
    
    auto value = isl_val_get_num_si(pointVal);
 
    container->push_back(value);
    return isl_stat_ok; 
  }    

isl_stat find_min_point_val(isl_point * pnt, void * user) {

  isl_val * pointVal = isl_point_get_coordinate_val(pnt, isl_dim_set, 0);
  
  auto value = isl_val_get_num_si(pointVal);
  
  int * currMin = (int*) user;
  
  if(*currMin > value) {
    *currMin = value;
  }
  return isl_stat_ok;
}
std::string lexico_constraints(int dim) {
    
        if(dim == 1) {
            std::string constraint = "i0 < o0";
            return constraint;
        }
       
        std::string lowerConstraints = lexico_constraints(dim-1);
        
        std::string newConstraint = "i" + std::to_string(dim-1) + " < o" + std::to_string(dim-1);

        for(int i = 0; i < dim; i++) {
            newConstraint += "and i" + std::to_string(i) + "= o" + std::to_string(i);
        }

        return lowerConstraints + " or " + newConstraint;

    }
    std::string get_array_name(SymAccess a) {
        
        static std::map<const void*,std::string> map;
        static int arrayCount = 0;
        
        
        auto find = map.find(a.view);
        if (find != map.end()) {
            
            return find->second;
        } else {
            std::stringstream newName;
            newName << "arr" << arrayCount;
            arrayCount++;
            auto newPair = std::pair<const void*,std::string>(a.view, newName.str());
            map.insert(newPair);
            return newName.str();
        }
    }

  std::string read_string(std::vector<std::vector<SymAccess>> accessLists, int dimensions) {
    std::string is = "";
    for(int i = 0; i < dimensions; i++) {
      is = is + "i" + std::to_string(i);
      if(i != dimensions-1) {
        is = is + ",";
      }
    }
    
    std::stringstream reads;
    reads << "{";
    for(long unsigned int i = 0; i < accessLists.size(); i++) {
      auto lst = accessLists.at(i);
      for(RAJA::SymAccess a : lst) {
        if(a.isRead) {
          std::string array_name = get_array_name(a);
          std::string loop_name = "L" + std::to_string(i);
          reads << "\t" << loop_name << "[" << is << "] -> " << array_name << "[" << a.access_string() << "];";
        }
      }
    }//accessLists.size()

    reads << "}";
    return reads.str();
  } 

  std::string write_string(std::vector<std::vector<SymAccess>> accessLists, int dimensions) {
    std::string is = ""; 
    for(int i = 0; i < dimensions; i++) {
      is = is + "i" + std::to_string(i);
      if(i != dimensions-1) {
        is = is + ",";
      }
    }
    
    std::stringstream writes;
    writes << "{";
    for(long unsigned int i = 0; i < accessLists.size(); i++) {
      auto lst = accessLists.at(i);
      for(RAJA::SymAccess a : lst) {
        if(a.isWrite) {
          std::string array_name = get_array_name(a);
          std::string loop_name = "L" + std::to_string(i);
          writes << "\t" << loop_name << "[" << is << "] -> " << array_name << "[" << a.access_string() << "];";
        }
      }
    }//accessLists.size()

    writes << "}";
    return writes.str();
  } 




    std::string read_string(std::vector<SymAccess> a1, std::vector<SymAccess> a2, int dimensions) {
        std::string bnds = "";
        std::string is = "";
        for(int i = 0; i < dimensions; i++) {
            is = is + "i" + std::to_string(i);
            bnds = bnds + "l" + std::to_string(i) + ",u" + std::to_string(i);
            if (i != dimensions-1){
               is  = is + ",";
               bnds += ",";
            }
        }
        
        std::stringstream reads;
         
        reads << /*"[" << bnds << "] ->*/"{";
        
        for(RAJA::SymAccess a : a1) {
            if(a.isRead) {
                std::string array_name = get_array_name(a);
                reads << "\tL1[" << is << "] -> " << array_name << "[" << a.access_string() << "];";
            }
            
        }
        for(RAJA::SymAccess a : a2) {
            if(a.isRead) {
                std::string array_name = get_array_name(a);
                reads << "\tL2[" << is << "] -> " << array_name << "[" << a.access_string() << "];";
            }
            
        }
        reads << "}";
        return reads.str();
    } 
   
    std::string write_string(std::vector<SymAccess> a1, std::vector<SymAccess> a2, int dimensions) {
        std::string bnds = "";
        std::string is = "";
        for(int i = 0; i < dimensions; i++) {
            is = is + "i" + std::to_string(i);
            bnds = bnds + "l" + std::to_string(i) + ",u" + std::to_string(i);
            if (i != dimensions-1){
               is  = is + ",";
               bnds += ",";
            }
        }
        
        std::stringstream writes;
         
        //writes << "[" << is << "] -> {";//"[" << bnds << "] -> {";
        writes << /*"[" << bnds << "] ->*/"{";
        
        for(RAJA::SymAccess a : a1) {
            if(a.isWrite) {
                std::string array_name = get_array_name(a);
                writes << "\tL1[" << is << "] -> " << array_name << "[" << a.access_string() << "];";
            }
            
        }
        for(RAJA::SymAccess a : a2) {
            if(a.isWrite) {
                std::string array_name = get_array_name(a);
                writes << "\tL2[" << is << "] -> " << array_name << "[" << a.access_string() << "];";
            }
            
        }
        writes << "}";
        return writes.str();
    } 

    std::string read_string_1d(std::vector<SymAccess> a1, std::vector<SymAccess> a2) {
        std::string lb = "l";
        std::string ub = "u";
        std::string bnds = lb + "," + ub;
        std::string i0 = "i0";
        
        std::stringstream reads;
        
        reads << "[" << bnds << "] -> {";
        
        for(RAJA::SymAccess a : a1) {
            if(a.isRead) {
                std::string array_name = get_array_name(a);
                reads << "\tL1[" << i0 << "] -> " << array_name << "[" << a.access_string() << "];";
            }
            
        }
        for(RAJA::SymAccess a : a2) {
            if(a.isRead) {
                std::string array_name = get_array_name(a);
                reads << "\tL2[" << i0 << "] -> " << array_name << "[" << a.access_string() << "];";
            }
            
        }
        
        reads << "}";
        return reads.str();
    }
    
    std::string read_string_2d(std::vector<SymAccess> a1, std::vector<SymAccess> a2) {
        std::string bnds = "l0,u0,l1,u1"; 
        std::string i0 = "i0";
        std::string i1 = "i1";
        std::string is = i0 + "," + i1;
        std::stringstream reads;
    

       reads << "[" << bnds << "] -> {";
      
       for(RAJA::SymAccess a : a1) {
           if(a.isRead) {
               std::string array_name = get_array_name(a);
               reads << "\tL1[" << is << "] -> " << array_name << "[" << a.access_string() << "];";
           }
       } 

       for(RAJA::SymAccess a : a2) {
           if(a.isRead) {
               std::string array_name = get_array_name(a);
               reads << "\tL2[" << is << "] -> " << array_name << "[" << a.access_string() << "];";
           }
       } 

       reads << "}";

       return reads.str();
    }
    std::string write_string_1d(std::vector<SymAccess> a1, std::vector<SymAccess> a2) {
        std::string lb = "l";
        std::string ub = "u";
        std::string bnds = lb + "," + ub;
        std::string i0 = "i0";
        
        std::stringstream writes;
        
        writes << "[" << bnds << "] -> {";
        
        for(RAJA::SymAccess a : a1) {
            if(a.isWrite ) {
                std::string array_name = get_array_name(a);
                writes << "\tL1[" << i0 << "] -> " << array_name << "[" << a.access_string() << "];";
            }
            
        }
        for(RAJA::SymAccess a : a2) {
            if(a.isWrite) {
                std::string array_name = get_array_name(a);
                writes << "\tL2[" << i0 << "] -> " << array_name << "[" << a.access_string() << "];";
            }
            
        }
        
        writes << "}";
        return writes.str();
    }
    int isl_can_fuse_2d(const char * readString, const char * writeString) {

        isl_ctx* ctx;
        isl_union_set* domain;
        isl_union_map* schedule, *read, *write, *schedule_inverse, *lex, *before;
        isl_union_map* raw, *war, *waw;
        
        ctx = isl_ctx_alloc();
        isl_printer* p = isl_printer_to_file(ctx, stdout);
        
        
        printf("\n Running isl_can_fuse_2d()\nRead String: %s\nWrite String: %s\n", readString,writeString);
        const char * domainString = "[l0,u0,l1,u1] -> { L1[i0,i1] : i0 < u0 and i0 >= l0 and i1 < u1 and i1 >= l1; L2[i0,i1] : i0 < u0 and i0 >= l0 and i1 < u1 and i1 >= l1;}";

        const char * scheduleString =  "[l0,u0,l1,u1] -> {L1[i0,i1] -> [0,i0,i1,0]; L2[i0,i1] -> [1,i0,i1,0];}";
        
        //const char * fusedScheduleString =  "[l0,u0,l1,u1] -> {L1[i0,i1] -> [0,i0,i1,0]; L2[i0,i1] -> [0,i0,i1,1];}";

        domain = isl_union_set_read_from_str(ctx,domainString);
        
        printf("\nDomain:\n");
        p = isl_printer_print_union_set(p, isl_union_set_copy(domain));
        
	//return 0;
        read = isl_union_map_read_from_str(ctx,readString);
        
        printf("\nReads, pre-intersect\n");
        p = isl_printer_print_union_map(p, read);
        
        read = isl_union_map_intersect_domain(read, isl_union_set_copy(domain));

        printf("\nRead, post-intersect:\n");
        p = isl_printer_print_union_map(p,read);

        write = isl_union_map_read_from_str(ctx, writeString);
        
        //printf("\nWrite, pre-intersect\n");
        //p = isl_printer_print_union_map(p,write);
        write = isl_union_map_intersect_domain(write, isl_union_set_copy(domain));
        
        //printf("\nWrite, post-intersection:\n");
        //p = isl_printer_print_union_map(p,write);
        
        schedule = isl_union_map_read_from_str(ctx, scheduleString);
        
        schedule_inverse = isl_union_map_reverse(isl_union_map_copy(schedule));
        //printf("\nSchedule:\n");
        //p = isl_printer_print_union_map(p,schedule);
        
        //printf("\nInverse Schedule:\n");
        //p = isl_printer_print_union_map(p,schedule_inverse);

        lex = isl_union_map_read_from_str(ctx,
                                          "{ [i0,i1,i2] -> [o0,o1,o2] : i0 < o0 or i0 = o0 and i1 < o1 or i0 = o0 and i1 = o1 and i2 < o2 }"
                                          );
        
        before = isl_union_map_apply_range(isl_union_map_copy(schedule), isl_union_map_copy(lex));
        before = isl_union_map_apply_range(before, isl_union_map_copy(schedule_inverse));
        
        //printf("\nBefore:\n");
        //p = isl_printer_print_union_map(p, before);
        
        isl_union_map * read_inverse = isl_union_map_reverse(isl_union_map_copy(read));
        isl_union_map * write_inverse = isl_union_map_reverse(isl_union_map_copy(write));
        
        raw = isl_union_map_apply_range(isl_union_map_copy(write), isl_union_map_copy(read_inverse));
        raw = isl_union_map_intersect(raw, isl_union_map_copy(before));
        
        waw = isl_union_map_apply_range(isl_union_map_copy(write), isl_union_map_copy(write_inverse));
        waw = isl_union_map_intersect(waw, isl_union_map_copy(before));
        
        war = isl_union_map_apply_range(isl_union_map_copy(read), isl_union_map_copy(write_inverse));
        war = isl_union_map_intersect(war, isl_union_map_copy(before));
        


        return 0;
    }
 
    int isl_can_fuse_1d(const char * readString, const char * writeString) {
        isl_ctx* ctx;
        isl_union_set* domain;
        isl_union_map* schedule, *read, *write, *schedule_inverse, *lex, *before;
        isl_union_map* raw, *war, *waw;
        
        ctx = isl_ctx_alloc();
        isl_printer* p = isl_printer_to_file(ctx, stdout);
        
        
        printf("\n Running isl_can_fuse()\nRead String: %s\nWrite String: %s\n", readString,writeString);
        const char * domainString = "[l,u] -> { L1[i] : i < u and i >= l; L2[i] : i < u and i >= l;}";
        const char * scheduleString =  "[l,u] -> {L1[i0] -> [0,i0,0]; L2[i0] -> [1,i0,0];}";
        
        const char * fusedScheduleString =  "[l,u] -> {L1[i0] -> [0,i0,0]; L2[i0] -> [0,i0,1];}";
        
        domain = isl_union_set_read_from_str(ctx,
                                             domainString
                                             );
        printf("\nDomain:\n");
        p = isl_printer_print_union_set(p, isl_union_set_copy(domain));
        
        read = isl_union_map_read_from_str(ctx, readString);
        
        
        printf("\nRead, pre-intersect:\n");
        p = isl_printer_print_union_map(p, read);
        read = isl_union_map_intersect_domain(read, isl_union_set_copy(domain));
        
        //printf("\nRead, post-intersect:\n");
        //p = isl_printer_print_union_map(p,read);
        
        write = isl_union_map_read_from_str(ctx,
                                            writeString
                                            //"[n] -> { S[k] -> C[k]; T[i, j] -> C[i + j];}"
                                            );
        
        //printf("\nWrite, pre-intersect\n");
        //p = isl_printer_print_union_map(p,write);
        write = isl_union_map_intersect_domain(write, isl_union_set_copy(domain));
        
        //printf("\nWrite, post-intersection:\n");
        //p = isl_printer_print_union_map(p,write);
        
        schedule = isl_union_map_read_from_str(ctx,
                                               //    "[n] -> { T[i, j] -> [1, i, j]; S[k] -> [0, k, 0];}"
                                               scheduleString);
        
        schedule_inverse = isl_union_map_reverse(isl_union_map_copy(schedule));
        //printf("\nSchedule:\n");
        //p = isl_printer_print_union_map(p,schedule);
        
        //printf("\nInverse Schedule:\n");
        //p = isl_printer_print_union_map(p,schedule_inverse);
        
        lex = isl_union_map_read_from_str(ctx,
                                          "{ [i0,i1,i2] -> [o0,o1,o2] : i0 < o0 or i0 = o0 and i1 < o1 or i0 = o0 and i1 = o1 and i2 < o2 }"
                                          );
        
        before = isl_union_map_apply_range(isl_union_map_copy(schedule), isl_union_map_copy(lex));
        before = isl_union_map_apply_range(before, isl_union_map_copy(schedule_inverse));
        
        //printf("\nBefore:\n");
        //p = isl_printer_print_union_map(p, before);
        
        isl_union_map * read_inverse = isl_union_map_reverse(isl_union_map_copy(read));
        isl_union_map * write_inverse = isl_union_map_reverse(isl_union_map_copy(write));
        
        raw = isl_union_map_apply_range(isl_union_map_copy(write), isl_union_map_copy(read_inverse));
        raw = isl_union_map_intersect(raw, isl_union_map_copy(before));
        
        waw = isl_union_map_apply_range(isl_union_map_copy(write), isl_union_map_copy(write_inverse));
        waw = isl_union_map_intersect(waw, isl_union_map_copy(before));
        
        war = isl_union_map_apply_range(isl_union_map_copy(read), isl_union_map_copy(write_inverse));
        war = isl_union_map_intersect(war, isl_union_map_copy(before));
        
        
        //printf("\nRAW:\n");
        //p = isl_printer_print_union_map(p,raw);
        //printf("\nWAW:\n");
        //p = isl_printer_print_union_map(p,waw);
        //printf("\nWAR:\n");
        //p = isl_printer_print_union_map(p,war);
        
        isl_union_map * fused_schedule = isl_union_map_read_from_str(ctx,
                                                                     fusedScheduleString);
        isl_union_map * fused_schedule_inverse = isl_union_map_reverse(isl_union_map_copy(fused_schedule));
        
        //printf("\nFused Schedule:\n");
        //p = isl_printer_print_union_map(p, fused_schedule);
        
        
        isl_union_map * fused_before = isl_union_map_apply_range(isl_union_map_copy(fused_schedule), isl_union_map_copy(lex));
        
        fused_before = isl_union_map_apply_range(fused_before, isl_union_map_copy(fused_schedule_inverse));
        
        //printf("\nFused Before\n";
        //p = isl_printer_print_union_map(p, fused_before);
        
        isl_bool respect1 = isl_union_map_is_subset(raw, isl_union_map_copy(fused_before));
        
        //printf("\nRespects RAW: %d\n", respect1);
        
        isl_bool respect2 = isl_union_map_is_subset(waw, isl_union_map_copy(fused_before));
        
        //printf("\nRespects WAW: %d\n", respect2);
        
        isl_bool respect3 = isl_union_map_is_subset(war, isl_union_map_copy(fused_before));
        
        //printf("\nRespects WAR: %d\n", respect3);
        
        int fusable = respect1 && respect2 && respect3;
        //printf("fusable: %d\n", fusable);
        
        return fusable;
        
    }
    
    int can_fuse_1d(std::vector<SymAccess> const & a1, std::vector<SymAccess> const & a2) {
      
        //std::cout << "can fuse\n"; 
        std::string read = read_string_1d(a1,a2);
        
        std::string write = write_string_1d(a1,a2);
        //std::cout << "readstring: " << read << "\nwritestring: " << write << "\n\n"; 
        int result =  isl_can_fuse_1d(read.c_str(), write.c_str());
        //std::cout << "can fuse: " << result << "\n";
        return result;
        return 0;
    }

    int can_fuse_2d(std::vector<SymAccess> const & a1, std::vector<SymAccess> const & a2) {
    
        std::string reads = read_string(a1,a2,2);
        std::string write = write_string(a1,a2,2);

        int result = isl_can_fuse_2d(reads.c_str(), write.c_str());

        return result;
    }



int amount_to_shift_ref(std::vector<SymAccess> const & a1, std::vector<SymAccess> const & a2) {
        isl_ctx* ctx;
        isl_union_set* domain;
        isl_union_map* schedule, *read, *write, *schedule_inverse, *lex, *before;
        isl_union_map* raw, *war, *waw;
        
        ctx = isl_ctx_alloc();
        isl_printer* p = isl_printer_to_file(ctx, stdout);
        
        std::cout << read_string_1d(a1,a2);
        std::cout << read_string_1d(a1,a2).c_str();
        auto rString = read_string_1d(a1,a2);
        auto wString = write_string_1d(a1,a2);

        std::cout << "\n\n" << rString;
        std::cout << rString.c_str();

        auto readString = rString.c_str();
        auto writeString = wString.c_str();

        printf("\n Running isl_can_fuse()\nRead String: %s\nWrite String: %s\n", readString,writeString);
        const char * domainString = "[l,u] -> { L1[i] : i < u and i >= l; L2[i] : i < u and i >= l;}";
        const char * scheduleString =  "[l,u] -> {L1[i0] -> [0,i0,0]; L2[i0] -> [1,i0,0];}";
        
        const char * fusedScheduleString =  "[l,u] -> {L1[i0] -> [0,i0,0]; L2[i0] -> [0,i0,1];}";
        
        domain = isl_union_set_read_from_str(ctx,
                                             domainString
                                             );
        printf("\nDomain:\n");
        p = isl_printer_print_union_set(p, isl_union_set_copy(domain));
        
        read = isl_union_map_read_from_str(ctx, readString);
        
        
        printf("\nRead, pre-intersect:\n");
        p = isl_printer_print_union_map(p, read);
        read = isl_union_map_intersect_domain(read, isl_union_set_copy(domain));
        
        //printf("\nRead, post-intersect:\n");
        //p = isl_printer_print_union_map(p,read);
        
        write = isl_union_map_read_from_str(ctx,
                                            writeString
                                            //"[n] -> { S[k] -> C[k]; T[i, j] -> C[i + j];}"
                                            );
        
        //printf("\nWrite, pre-intersect\n");
        //p = isl_printer_print_union_map(p,write);
        write = isl_union_map_intersect_domain(write, isl_union_set_copy(domain));
        
        //printf("\nWrite, post-intersection:\n");
        //p = isl_printer_print_union_map(p,write);
        
        schedule = isl_union_map_read_from_str(ctx,
                                               //    "[n] -> { T[i, j] -> [1, i, j]; S[k] -> [0, k, 0];}"
                                               scheduleString);
        
        schedule_inverse = isl_union_map_reverse(isl_union_map_copy(schedule));
        //printf("\nSchedule:\n");
        //p = isl_printer_print_union_map(p,schedule);
        
        //printf("\nInverse Schedule:\n");
        //p = isl_printer_print_union_map(p,schedule_inverse);
        
        lex = isl_union_map_read_from_str(ctx,
                                          "{ [i0,i1,i2] -> [o0,o1,o2] : i0 < o0 or i0 = o0 and i1 < o1 or i0 = o0 and i1 = o1 and i2 < o2 }"
                                          );
        
        before = isl_union_map_apply_range(isl_union_map_copy(schedule), isl_union_map_copy(lex));
        before = isl_union_map_apply_range(before, isl_union_map_copy(schedule_inverse));
        
        //printf("\nBefore:\n");
        //p = isl_printer_print_union_map(p, before);
        
        isl_union_map * read_inverse = isl_union_map_reverse(isl_union_map_copy(read));
        isl_union_map * write_inverse = isl_union_map_reverse(isl_union_map_copy(write));
        
        raw = isl_union_map_apply_range(isl_union_map_copy(write), isl_union_map_copy(read_inverse));
        raw = isl_union_map_intersect(raw, isl_union_map_copy(before));
        
        waw = isl_union_map_apply_range(isl_union_map_copy(write), isl_union_map_copy(write_inverse));
        waw = isl_union_map_intersect(waw, isl_union_map_copy(before));
        
        war = isl_union_map_apply_range(isl_union_map_copy(read), isl_union_map_copy(write_inverse));
        war = isl_union_map_intersect(war, isl_union_map_copy(before));
        
        isl_union_map * fused_schedule = isl_union_map_read_from_str(ctx,
                                                                     fusedScheduleString);
        isl_union_map * fused_schedule_inverse = isl_union_map_reverse(isl_union_map_copy(fused_schedule));
        
        //printf("\nFused Schedule:\n");
        //p = isl_printer_print_union_map(p, fused_schedule);
        
        
        isl_union_map * fused_before = isl_union_map_apply_range(isl_union_map_copy(fused_schedule), isl_union_map_copy(lex));
        
        fused_before = isl_union_map_apply_range(fused_before, isl_union_map_copy(fused_schedule_inverse));

        isl_union_map * raw_remaining = isl_union_map_subtract(raw, isl_union_map_copy(fused_before));
        isl_union_map * waw_remaining = isl_union_map_subtract(waw, isl_union_map_copy(fused_before));
        isl_union_map * war_remaining = isl_union_map_subtract(war, isl_union_map_copy(fused_before));

        printf("\nraw remaining\n");
        p = isl_printer_print_union_map(p, raw_remaining);
        printf("\nwaw remaining\n");
        p = isl_printer_print_union_map(p, waw_remaining);
        printf("\nwar remaining\n");
        p = isl_printer_print_union_map(p, war_remaining);
       
        return 0;
    }
} // namespace RAJA
