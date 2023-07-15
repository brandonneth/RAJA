#ifndef RAJA_FormatDecisions_HPP
#define RAJA_FormatDecisions_HPP

#include <typeinfo>
#include <cmath>
#include <chrono>
#include "RAJA/loopchain/FormatDecisions/Model/Constraint.hpp"
#include "RAJA/loopchain/FormatDecisions/Model/DecisionVariable.hpp"
#include "RAJA/loopchain/FormatDecisions/Model/ISLTranslation.hpp"
#include "RAJA/loopchain/FormatDecisions/Conversion.hpp"
#include "RAJA/loopchain/FormatDecisions/Enumerate.hpp"
#include "RAJA/loopchain/Utils.hpp"
namespace RAJA
{


template <typename T>
T get_array_name(SymAccess a);


template <typename KernelTuple, typename... ViewRefs>
struct FormatDecisions {
  constexpr static idx_t NumViews = sizeof...(ViewRefs);
  constexpr static idx_t NumKernels = camp::tuple_size<KernelTuple>::value;
  constexpr static idx_t NumMoments = NumKernels + 2;

  using LayoutTupleType =
      tuple<typename std::remove_reference<ViewRefs>::type::layout_type...>;

  tuple<ViewRefs...> views;
  KernelTuple knls;

  LayoutTupleType input_formats;

  std::array<LayoutTupleType, NumMoments> moments;
  std::array<std::array<bool, NumViews>, NumMoments> user_decided; 

  idx_t setup_time = 0;
  idx_t solve_time = 0;
  idx_t isl_time = 0;
  idx_t space_time = 0;
  idx_t map_time = 0;
  idx_t cost = 0;

  template <idx_t...Is>
  FormatDecisions(tuple<ViewRefs...> _views, KernelTuple _knls, camp::idx_seq<Is...>) : 
    views(_views), 
    knls(_knls) {
    
    moments[0] = make_tuple(get<Is>(views).get_layout()...);

    for(int i0 = 0; i0 < NumMoments; i0++ ) {
      for(auto i1 = 0; i1 < NumViews; i1++) {
        user_decided[i0][i1] = false;
      }
      moments[i0] = make_tuple(get<Is>(views).get_layout()...);

    }
    for(int i = 0; i < NumViews; i++) {
      user_decided[0][i] = true;
    }
  }

  template <camp::idx_t I, typename ViewType, typename InternalViewType>
  camp::idx_t view_index_helper(ViewType & view, InternalViewType & internal_view) {
    if constexpr (std::is_same<ViewType, InternalViewType>::value) {
      if (&view == &internal_view) {
        return I;
      } else {
        return -1;
      }
    } else {
      return -1;
    }
  }
  template <typename ViewType>
  camp::idx_t view_index(ViewType & view) {
    auto helper = [&]
    <camp::idx_t...Is>
    (ViewType & view, camp::idx_seq<Is...>) {
      return vmax(view_index_helper<Is>(view, get<Is>(views))...);
    };
    return helper(view, idx_seq_for(views));
  }

  template <typename KernelType>
  camp::idx_t kernel_index(KernelType) {
    auto helper1 = [&]
    <typename Knl>
    (Knl, int i) {
      if constexpr (std::is_same<Knl,KernelType>::value) {
         return i;
      } else {
         return -1;
      }
    };

    auto helper2 = [&]
    <camp::idx_t...Is>
    (camp::idx_seq<Is...>) {
      return vmax(helper1(camp::get<Is>(knls), Is)...);
    };

    return helper2(idx_seq_for(knls));
  }
  
  template <typename ViewType, typename LayoutType>
  int update_layout(LayoutTupleType & layouts, ViewType & view, LayoutType layout) {
    auto helper1 = [&]
    <camp::idx_t I>
    () {
      if constexpr (std::is_same<ViewType&,decltype(get<I>(this->views))>::value) {
        if(view_index(view) == I) {
          get<I>(layouts) = layout;
          return 1;
        }
      }
      return 0;
    };

    auto helper2 = [&]
    <camp::idx_t...Is>
    (camp::idx_seq<Is...>) {
      return (helper1.template operator()<Is>() + ...);
    };

    auto success = helper2(idx_seq_for(views));
    return success;
  }

  template <idx_t MomentIndex, typename ViewType>
  auto get_layout_choice(ViewType & view) {
    auto layouts = moments[MomentIndex];

    std::vector<idx_t> perm;

    auto helper1 = [&]
    <camp::idx_t I>
    () {
      if constexpr (std::is_same<ViewType&,decltype(get<I>(this->views))>::value) {
        if(view_index(view) == I) {
          auto p = layout_to_perm(get<I>(layouts));
          for(unsigned long int i = 0; i < p.size(); i++) {
            perm.push_back(p[i]);
          }
        }
      }
      return 0;
    };

    auto helper2 = [&]
    <camp::idx_t...Is>
    (camp::idx_seq<Is...>) {
      return (helper1.template operator()<Is>() + ...);
    };

    helper2(idx_seq_for(views));
    
    return perm;
  }

  template <idx_t CurrIdx, typename ViewType>
  auto user_constraints(ViewType & view) {
    using LayoutType = typename ViewType::layout_type;
    constexpr camp::idx_t NumDims = LayoutType().n_dims;
    
    if constexpr (CurrIdx >= NumMoments) {
      return new EqualityConstraint<NumDims>(new Constant<NumDims>(1), new Constant<NumDims>(1));
    } else {
      auto next_constraint = user_constraints<CurrIdx+1>(view);
    
      if (!user_decided[CurrIdx][view_index(view)]) {
        auto this_constraint = new EqualityConstraint<NumDims>(new Constant<NumDims>(1), new Constant<NumDims>(1));
        return new AndConstraint<NumDims>(next_constraint, this_constraint); 
      } else {
        auto chosen_perm = get_layout_choice<CurrIdx>(view);
        auto chosen_format_variable = new FormatDecision<NumDims>(view_name(view), CurrIdx, chosen_perm);
        auto this_constraint = new EqualityConstraint<NumDims>(new Constant<NumDims>(1), chosen_format_variable);
        return new AndConstraint<NumDims>(next_constraint, this_constraint);
      }
    }
  }

  template <typename ViewType>
  auto view_name(ViewType & view) {
    auto idx = view_index(view);
    std::stringstream s;
    if (idx == -1) {
      s << "false_array";
    } else {
      s << "arr" << idx;
    }
    return s.str();
  }

  template <idx_t NumDims>
  idx_t estimate_format_cost(std::array<idx_t, NumDims> layout_order) {
    static std::map<std::array<idx_t,NumDims>,idx_t> cache;

    if( cache[layout_order] != 0) {
      return cache[layout_order];
    }
    //std::cerr << "no cache\n";
    using V =  View<double, Layout<NumDims>>;
    constexpr camp::idx_t n = 128;
    camp::idx_t s = 1;
    std::array<camp::idx_t, NumDims> sizes;
    for(auto i = 0; i < NumDims; i++) {
      sizes[i] = n;
      s *= n;
    }

    Layout<NumDims> in_layout = make_permuted_layout(sizes, layout_order);

    auto data1 = new double[s];
    auto data2 = new double[s];
    V in1(data1, in_layout);
    V in2(data2, in_layout);

    std::srand(0);
    for(int i = 0; i < s; i++) {
      data1[i] = (double) std::rand();
      data2[i] = (double) std::rand();
    }
    using kpol = typename nest<NumDims>::pol;

    
    auto l = f_lambda(in1, in2);
    auto seg = RangeSegment(0,n);
    auto segs = tuple_repeat<NumDims>(seg);
    auto knl = make_kernel<kpol>(segs, l);

    auto start =std::chrono::high_resolution_clock::now();
    knl();
    auto stop =std::chrono::high_resolution_clock::now();

    auto t =std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start).count();
    cache[layout_order] = t;
    //std::cerr << "compcost: " << cache[layout_order] << "\n";
    return cache[layout_order];
  }

  template <idx_t KnlIdx, typename ViewType>
  auto estimate_format_costs(ViewType & view) {
    using LayoutType = typename ViewType::layout_type;
    constexpr camp::idx_t NumDims = LayoutType().n_dims;

    if constexpr (KnlIdx >= NumKernels) {
      return new SumExpression<NumDims>();
    } else {
      SumExpression<NumDims> * term_sum = estimate_format_costs<KnlIdx+1>(view);
      
      auto knl = get<KnlIdx>(knls);
      using KernelType = decltype(knl);    
      using KernelPolicy = typename KernelType::KPol;
      auto accesses = knl.execute_symbolically();
    
      for(auto access: accesses) {
        if(access.view == view.get_data()) {
          auto all_perms = enumerate_perms<NumDims>();
          for (auto perm : all_perms) {
            auto format_var = new FormatDecision<NumDims>(view_name(view), KnlIdx+1, perm);
          
            auto normalized_layout_perm = access.template normalized_layout_perm<KernelPolicy,NumDims>(perm);
            //std::cerr << "estimating cost of " << format_var->rep() << "\n";
            auto cost = estimate_format_cost<NumDims>(normalized_layout_perm);
            unsigned long coef = cost * knl.num_iterations();
            auto cost_term = new MultExpression<NumDims>(format_var, new Constant<NumDims>(coef));
            term_sum->add(cost_term);
          }
        }
      } 
      return term_sum;
    }
  }

  template <typename ViewType>
  auto finalize_view(ViewType & view) {
    auto setup_start = std::chrono::high_resolution_clock::now();
    using LayoutType = typename ViewType::layout_type;
    constexpr camp::idx_t NumDims = LayoutType().n_dims;
    
    SumExpression<NumDims> * cost_function = estimate_format_costs<0>(view);

    auto decision_variables = all_decision_variables<NumDims>(view_name(view), NumKernels);
    auto base_constraints = all_constraints<NumDims>(view_name(view), NumKernels); 
    auto extra_constraints = user_constraints<0>(view);
    auto all_constraints = new AndConstraint<NumDims>(base_constraints, extra_constraints);
    //TODO: this is a hack. The estimate_cost should probably take a policy or something and then that should be uniformized instead of doing format costs here. 
    for(auto v : decision_variables->variables) {
      auto estimate = v->estimate_cost() / 1000;
      unsigned long coef = estimate * view.size();
      auto c = new Constant<NumDims>(coef);
      MultExpression<NumDims> * cost_term = new MultExpression<NumDims>(v, c);
      cost_function->add(cost_term);
    }
    
    auto setup_stop = std::chrono::high_resolution_clock::now();
    setup_time += std::chrono::duration_cast<std::chrono::microseconds>(setup_stop - setup_start).count();

    auto isl_start = std::chrono::high_resolution_clock::now();
    
    //std::cerr << "cost function: " << cost_function->rep() << "\n";
    //std::cerr << "Setting up and solving ILP problem\n";
    //std::cerr << "User Constraints: " << extra_constraints->rep() << "\n";
    //auto start = std::chrono::system_clock::now();
    isl_ctx * ctx = isl_ctx_alloc();
     

    auto space_start = std::chrono::system_clock::now();
    isl_union_set * all_choices = decision_space(ctx, decision_variables, all_constraints);
    auto space_end = std::chrono::system_clock::now();
    space_time += std::chrono::duration_cast<std::chrono::microseconds>(space_end - space_start).count();
    //auto time_space = std::chrono::system_clock::now();

    auto map_start = std::chrono::system_clock::now();
    isl_union_map * objective_map = objective_function(ctx, decision_variables, all_constraints, cost_function);
    auto map_end = std::chrono::system_clock::now();
    map_time += std::chrono::duration_cast<std::chrono::microseconds>(map_end - map_start).count();
 
    auto solve_start = std::chrono::high_resolution_clock::now();
    isl_union_set * applied = isl_union_set_apply(all_choices, objective_map);
    isl_union_set * best_score = isl_union_set_lexmin(applied);
    isl_point * choice = isl_union_set_sample_point(best_score);
    auto solve_stop = std::chrono::high_resolution_clock::now();
    solve_time += std::chrono::duration_cast<std::chrono::microseconds>(solve_stop - solve_start).count();
/*
    std::cerr << "Made choice\n";
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> solve_time = end-start;
    std::chrono::duration<double> apply_time = time_applied-start;
    std::chrono::duration<double> create_time = time_create-start;
    std::chrono::duration<double> space_time = time_space-start;
    std::chrono::duration<double> alloc_time = time_space-start;
    std::cerr << "Elasped Time: " << solve_time.count() << "\n";
    std::cerr << "Time to create the ctx: " << alloc_time.count() << "\n";
    std::cerr << "Time to create the set : " << space_time.count() << "\n";
    std::cerr << "Time to create the set and  map: " << create_time.count() << "\n";
    std::cerr << "Time to create and apply the map: " << apply_time.count() << "\n";*/
    if(isl_point_is_void(choice)) {
      std::cerr << "No point found. Not making any modeling choices.\n";
     // isl_ctx_free(ctx); 
      return false;
    } 

    isl_val * score = isl_point_get_coordinate_val(choice, isl_dim_out, 0);
    cost += isl_val_get_num_si(score);
    int num_format_decisions = NumMoments * tgamma(NumDims+1);
    for(int i = 0; i < num_format_decisions; i++) {
      Term<NumDims> * var = decision_variables->variables.at(i);
      //std::cout << "Getting output value for " << var->rep() << ": ";
      isl_val * val_ = isl_point_get_coordinate_val(choice, isl_dim_out, i+1);
      long val = isl_val_get_num_si(val_);
       //std::cout << val << "\n";

      if(val == 1) {
        FormatDecision<NumDims> * vv = dynamic_cast<FormatDecision<NumDims>*>(var);
        auto perm = vv->permutation;
        auto idx = vv->idx;
        std::array<idx_t, NumDims> sizes;
        for(auto s = 0; s < NumDims; s++) {
          sizes[s] = view.get_layout().sizes[s];
        }
        auto newlayout = make_permuted_layout(sizes, perm);
        update_layout(moments[idx], view, newlayout);

      }
    }

    //isl_ctx_free(ctx); 
    return true;
  }

  template <typename ViewType, typename LayoutType, typename KernelType>
  bool set_format_for(ViewType & view, LayoutType layout, KernelType knl) {
    camp::idx_t v_idx = view_index(view);
    camp::idx_t k_idx = kernel_index(knl);
    if(v_idx == -1) {
      std::cerr << "View not found in tuple of views for this decision object\n";
      return false;
    }
    if(k_idx == -1) {
      std::cerr << "Kernel not found in this decision object\n";
      return false;
    }
    auto success = update_layout(moments[k_idx+1], view, layout);
    user_decided[k_idx+1][v_idx] = true;
    return success;

  }

  template <typename ViewType, typename LayoutType>
  bool set_output_format(ViewType & view, LayoutType layout) {
    camp::idx_t v_idx = view_index(view);
    if(v_idx == -1) {
      std::cerr << "View not found in tuple of views for this decision object\n";
      return false;
    }
    auto success = update_layout(moments[NumKernels+1], view, layout);
    user_decided[NumKernels+1][v_idx] = true;
    return success;
  }

  auto finalize() {
    isl_time = 0;
    setup_time = 0;
    solve_time = 0;
    map_time = 0;
    space_time = 0;
    auto view_finalizer = [&]
    <camp::idx_t...Is>
    (camp::idx_seq<Is...>) {
      camp::sink(finalize_view(get<Is>(views))...);
    };
    view_finalizer(idx_seq_for(views));

    auto moment_generator = [&]
    <camp::idx_t...Is>
    (camp::idx_seq<Is...>) {
      return make_tuple(permute_views(views, moments[Is])...);
    };
    auto conversions = moment_generator(idx_seq_from_to<1,NumMoments>());

    auto stitch = [&]
    <camp::idx_t CurrStitch>
    () {
      return tuple_cat(get<CurrStitch>(conversions), make_tuple(get<CurrStitch>(knls)));
    };
    auto stitcher = [&]
    <camp::idx_t...Is>
    (camp::idx_seq<Is...>) {
      auto final_conversion = get<NumKernels>(conversions);
      auto all_but_final = tuple_cat(stitch.template operator()<Is>()...);
      return tuple_cat(all_but_final, final_conversion);
    };
    auto knls_and_conversions= stitcher(idx_seq_for(knls));

    return loop_chain(knls_and_conversions);
  }

  void lock() {
    for(auto i= 0; i < NumMoments; i++) {
      for(auto j = 0; j < NumViews; j++) {
        user_decided[i][j] = true;
      }
    }
  }

  template <idx_t I>
  void print_to_stream_helper(std::ostream & o, LayoutTupleType t) {
    if constexpr (I >= NumViews) {
      return;
    } else {
      auto curr_layout = get<I>(t);

      auto perm = layout_to_perm(curr_layout);
      auto s = "";
      o << "(";
      for(int i = 0; i < perm.size(); i++) {
        o << s << perm[i];
        s  = ",";
      }
        o << ")";
      print_to_stream_helper<I+1>(o,t);
    }
  }
  auto print_to_stream(std::ostream & o) {
    o << "[";
    for(auto i = 0; i < NumMoments; i++) {
      print_to_stream_helper<0>(o, moments[i]);
    }
    o << "]";
  }
  
  template <typename KnlTuple, camp::idx_t I, camp::idx_t...Is>
  camp::tuple<idx_t, idx_t> te_helper(KnlTuple knlTuple, camp::idx_seq<I, Is...>, camp::tuple<idx_t, idx_t> totals) {
    

    auto start = std::chrono::high_resolution_clock::now();
    camp::get<I>(knlTuple)();
    auto stop = std::chrono::high_resolution_clock::now();
    camp::idx_t elapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

    
    if constexpr ((I+1) % (NumViews+1) != 0) {
      camp::get<0>(totals) += elapsed;
    } else {
      camp::get<1>(totals) += elapsed;
    } 
    if constexpr (sizeof...(Is) == 0) {
      return totals;
    } else {
      return te_helper(knlTuple, camp::idx_seq<Is...>{}, totals);
    }
  }
  camp::tuple<idx_t, idx_t> time_execution() {
    auto loopchain = finalize();
    auto knlTuple = loopchain.knlTuple;

    return te_helper(knlTuple, idx_seq_for(knlTuple), camp::tuple<idx_t, idx_t>{0,0});
  }
};

template <typename ViewRefTuple, typename...KernelTypes>
auto format_decisions(ViewRefTuple views, KernelTypes...knls) {
  using KernelTuple = decltype(make_tuple(knls...));

  auto helper = [=]
  <typename...ViewRefs, camp::idx_t...Is>
  (tuple<ViewRefs...> views) {
    return FormatDecisions<KernelTuple,ViewRefs...>(views, make_tuple(knls...), idx_seq_for(views));
  };

  return helper(views);
}

}  // namespace RAJA

#endif
