array[,] vector assemble_get_caught_and_reported_step(
  array[,] vector selected_weighted_fishing_step,
  vector reporting_rate
) {
  // Get dimensions
  int T = dims(selected_weighted_fishing_step)[1];
  int L = dims(selected_weighted_fishing_step)[2];
  int S = dims(selected_weighted_fishing_step)[3];
  // Declare values
  array[T, L] vector[S] get_caught_and_reported_step;
  // Populate
  for (t in 1:T) {
    for (l in 1:L) {
      get_caught_and_reported_step[t, l] = reporting_rate
      .* (1.0 - exp(-selected_weighted_fishing_step[t, l]));
    }
  }
  // Return value
  return get_caught_and_reported_step;
}

array[] vector assemble_selectivity (
  array[] vector selectivity_short,
  int L,
  int S
  ) {
  // Declare values
  array[L] vector[S] selectivity;
  // Populate values
  for (l in 1:L) {
    if (l == L) {
      selectivity[l] = rep_vector(1.0, S);
    } else {
      selectivity[l] = selectivity_short[l];
    }
  }
  // Return value
  return selectivity;
}

array[,] vector assemble_selected_weighted_fishing_step (
  array[] vector fishing_rate,
  // array[] vector fishing_weight,
  array[] vector selectivity,
  int T,
  int K
) {
  // Get dimensions
  int Years = dims(fishing_rate)[1];
  int L = dims(selectivity)[1];
  int S = dims(selectivity)[2];
  // Initiate values
  array[T, L] vector[S] selected_weighted_fishing_step;
  real years_per_step = 1 / (1.0 * K);
  int t = 1;
  // Populate value
  for (year in 1:Years) {
    for (k in 1:K) {
      for (l in 1:L) {
        selected_weighted_fishing_step[t, l] = selectivity[l]
        * years_per_step
        //.* fishing_weight[k]
        .* fishing_rate[year];
      }
      // Increment index
      t += 1;
    }
  }
  // Return value
  return selected_weighted_fishing_step;
}

array[] int assemble_simplex_dims (array[,] int mindex) {
  // Get dimensions
  int S = dims(mindex)[1];
  int A = 6;
  // Declare values
  array[A] int simplex_dimensions = rep_array(0, A);
  int row_x_sum;
  // Populate simplex dimensions
  for (s0 in 1:S) {
    row_x_sum = sum(mindex[s0]);
    if (row_x_sum > 0) {
      if (row_x_sum > A) {
        reject("row_x_sum: ", row_x_sum);
      }
      simplex_dimensions[row_x_sum] += 1;
    }
  }
  // Return value
  return simplex_dimensions;
}

array[,] matrix assemble_survive_then_move_step (
  array[] matrix movement_step,
  array[,] vector selected_weighted_fishing_step,
  vector natural_mortality_plus_loss_step
){
  // Get dimensions
  int T = dims(selected_weighted_fishing_step)[1];
  int L = dims(selected_weighted_fishing_step)[2];
  int S = dims(selected_weighted_fishing_step)[3];
  // Declare values
  array[T, L] matrix[S, S] survive_then_move_step;
  // Populate
  for (t in 1:T) {
    for (l in 1:L) {
      survive_then_move_step[t, l] = diag_pre_multiply(
        // Survival step
        exp(
          -selected_weighted_fishing_step[t, l]
          - natural_mortality_plus_loss_step
        ),
        // Movement step
        movement_step[l]
      );
    }
  }
  // Return value
  return survive_then_move_step;
}

array[,] vector assemble_tags_released (array[,,,,] int tags) {
  // Get dimensions
  int T = dims(tags)[1]; // Here T = model T - 1
  int L = dims(tags)[3];
  int S = dims(tags)[4];
  // Declare values
  array[T, L] vector[S] tags_released;
  // Populate tags released
  for (t in 1:T) { // Here T = model T - 1
    for (l in 1:L) { // Released size
      for (s0 in 1:S) { // Released region
        tags_released[t, l, s0] = tags[t, 1, l, s0, s0] * 1.0;
      }
    }
  }
  // Return tags released
  return tags_released;
}

array [,,,,] int assemble_tags_transpose (array[,,,,] int tags) {
  // Get dimensions
  int T = dims(tags)[1]; // Here T = model T - 1
  int D = dims(tags)[2];
  int L = dims(tags)[3];
  int S = dims(tags)[4];
  // Declare values
  array[T, D, L, S, S] int tags_transpose;
  // Populate tag array
  for (t in 1:T) { // Here T = model T - 1
    for (d in 1:D) {
      for (l in 1:L) {
        for (s0 in 1:S) {
          for (s in 1:S) {
            tags_transpose[t, d, l, s, s0] = tags[t, d, l, s0, s];
          }
        }
      }
    }
  }
  // Return tag array
  return tags_transpose;
}

/**
* Assemble Matrices That Indicate Possible Movement at Each Duration Step
*
* @param movement_matrix, a matrix of dimension [S, S]
* @param D, an integer giving the maximum duration at large in steps
*
* @return an array of dimension [D] holding matrices of dimension [S, S]
*
* The argument mindex is a square integer array of zeros and ones indicating
* that movement is permitted (one) or not permitted (zero) between a given
* source region (row) and destination region (column) in one model step.
*/
array[] matrix assemble_movement_possible (
  matrix mmatrix,
  int D
) {
  // Get dimensions
  int S = dims(mmatrix)[1];
  // Initialize values
  array[D] matrix[S, S] movement_possible;
  // Populate movement possible
  movement_possible[1] = rep_matrix(0.0, S, S);
  movement_possible[2] = mmatrix;
  // Iterate higher indexes
  for (d in 3:D) {
    if (min(movement_possible[d - 1]) > 0.0) {
      movement_possible[d] = movement_possible[d - 1];
    } else {
      movement_possible[d] = movement_possible[d - 1] * mmatrix;
    }
  }
  // Return movement possible
  return movement_possible;
}

array[] matrix assemble_movement_step (
  array[,] vector m1, // [L, ]
  array[,] vector m2, // [L, ]
  array[,] vector m3, // [L, ]
  array[,] vector m4, // [L, ]
  array[,] vector m5, // [L, ]
  array[,] vector m6, // [L, ]
  array[,] int mindex, // [S, S]
  int L
) {
  // Get dimensions
  int S = dims(mindex)[1];
  int A = 6;
  // Declare values
  array[L] matrix[S, S] movement_step = rep_array(rep_matrix(0.0, S, S), L);
  array[L, A] int index = rep_array(0, L, A); // Simplex array row index
  int row_x_sum; // Movement index row sum
  int column; // Simplex index (column)
  // Populate movement step
  for (l in 1:L) {
    for (s0 in 1:S) {
      row_x_sum = sum(mindex[s0]);
      if (row_x_sum > 0) {
        index[l, row_x_sum] += 1;
        column = 0;
        for (s in 1:S) {
          if (mindex[s0, s] == 1) {
            column += 1;
            if (row_x_sum == 1) {
              movement_step[l, s0, s] = m1[l, index[l, row_x_sum], column];
            } else if (row_x_sum == 2) {
              movement_step[l, s0, s] = m2[l, index[l, row_x_sum], column];
            } else if (row_x_sum == 3) {
              movement_step[l, s0, s] = m3[l, index[l, row_x_sum], column];
            } else if (row_x_sum == 4) {
              movement_step[l, s0, s] = m4[l, index[l, row_x_sum], column];
            } else if (row_x_sum == 5) {
              movement_step[l, s0, s] = m5[l, index[l, row_x_sum], column];
            } else if (row_x_sum == 6) {
              movement_step[l, s0, s] = m6[l, index[l, row_x_sum], column];
            } else {
              reject("row_x_sum: ", row_x_sum);
            }
          }
        }
      }
    }
  }
  // Return movement step
  return movement_step;
}

// INFO: Updated above 2024-08-23
// array[,] matrix assemble_survive_then_move_step (
//   array[] matrix movement_step,
//   array[] vector fishing_rate,
//   // array[] vector fishing_weight,
//   array[] vector selectivity,
//   vector natural_mortality_rate,
//   real ongoing_loss_rate,
//   int T,
//   int K
// ) {
//   // Get dimensions
//   int Years = dims(fishing_rate)[1];
//   int L = dims(selectivity)[1];
//   int S = dims(selectivity)[2];
//   // Initiate values
//   array[T, L] matrix[S, S] survive_then_move_step;
//   array[T, L] vector[S] survival_step;
//   real years_per_step = 1 / (1.0 * K);
//   int t = 1;
//   // Populate survive then move step
//   for (year in 1:Years) {
//     for (k in 1:K) {
//       for (l in 1:L) {
//         survive_then_move_step[t, l] = diag_pre_multiply(
//           // Survival step
//           exp(
//             // -selectivity[l] .* fishing_weight[k] .* fishing_rate[year]
//             -selectivity[l] .* fishing_rate[year] * years_per_step
//             - years_per_step * (natural_mortality_rate + ongoing_loss_rate)
//           ),
//           // Movement step
//           movement_step[l]
//         );
//       }
//       // Increment index
//       t += 1;
//     }
//   }
//   // Return value
//   return survive_then_move_step;
// }

// INFO: Updated above 2024-08-23
// array[,] vector assemble_survival_step (
//   array[] vector fishing_rate,
//   // array[] vector fishing_weight,
//   array[] vector selectivity,
//   vector natural_mortality_rate,
//   real ongoing_loss_rate,
//   int T,
//   int K
// ) {
//   // Get dimensions
//   int Years = dims(fishing_rate)[1];
//   int L = dims(selectivity)[1];
//   int S = dims(selectivity)[2];
//   // Initiate values
//   array[T, L] vector[S] survival_step;
//   real years_per_step = 1 / (1.0 * K);
//   int t = 1;
//   // Populate survival step
//   for (year in 1:Years) {
//     for (k in 1:K) {
//       for (l in 1:L) {
//         survival_step[t, l] = exp(
//           // -selectivity[l] .* fishing_weight[k] .* fishing_rate[year]
//           -selectivity[l] .* fishing_rate[year] * years_per_step
//           - years_per_step * (natural_mortality_rate + ongoing_loss_rate)
//         );
//       }
//       // Increment index
//       t += 1;
//     }
//   }
//   // Return value
//   return survival_step;
// }

// INFO: Updated above 2024-08-23
// array[,,] vector assemble_survival_step (
//   array[] vector fishing_step,
//   // array[] vector fishing_weight,
//   array[] vector selectivity,
//   vector natural_mortality_step,
//   real ongoing_loss_step,
//   int K,
//   int L
// ) {
//   // Get dimensions
//   int Years = dims(fishing_step)[1];
//   int S = dims(fishing_step)[2];
//   // Initialize values
//   array[Years, K, L] vector[S] survival_step;
//   // Populate survival step
//   for (year in 1:Years) {
//     for (k in 1:K) {
//       for (l in 1:L) {
//         survival_step[year, k, l] = exp(
//           -fishing_step[year] .* selectivity[l] // .* fishing_weight[k] .* selectivity[l]
//           - natural_mortality_step
//           - ongoing_loss_step
//         );
//       }
//     }
//   }
//   // Return survival step
//   return survival_step;
// }

// INFO: Updated above 2024-08-23
// array[,] matrix assemble_transition_step (
//   array[] matrix movement_step,
//   array[,,] vector survival_step
// ) {
//   // Get dimensions
//   int Years = dims(survival_step)[1];
//   int K = dims(survival_step)[2];
//   int L = dims(survival_step)[3];
//   int S = dims(survival_step)[4];
//   int T = Years * K;
//   // Declare values
//   array[T, L] matrix[S, S] transition_step;
//   int t = 1;
//   // Populate transition step
//   for (year in 1:Years) {
//     for (k in 1:K) {
//       for (l in 1:L) {
//         transition_step[t, l] = diag_pre_multiply( // A_n = A_{t-1}S_{t-1}\Gamma
//           survival_step[year, k, l],
//           movement_step[l]
//         );
//       }
//       t += 1;
//     }
//   }
//   // Return transition_step
//   return transition_step;
// }

// INFO: Updated above 2024-08-23
// array[,] vector assemble_observation_step (
//   array[] vector fishing_rate,
//   // array[] vector fishing_weight,
//   array[] vector selectivity,
//   vector reporting_rate,
//   int T,
//   int K
// ) {
//   // Get dimensions
//   int Years = dims(fishing_rate)[1];
//   int L = dims(selectivity)[1];
//   int S = dims(selectivity)[2];
//   // Initiate values
//   array[T, L] vector[S] observation_step;
//   real years_per_step = 1 / (1.0 * K);
//   int t = 1;
//   // Populate survival step
//   for (year in 1:Years) {
//     for (k in 1:K) {
//       for (l in 1:L) {
//         observation_step[t, l] = reporting_rate
//         // .* (1.0 - exp(-selectivity[l] .* fishing_weight[k] .* fishing_rate[year]))
//         .* (1.0 - exp(-selectivity[l] .* fishing_rate[year] * years_per_step));
//       }
//       // Increment index
//       t += 1;
//     }
//   }
//   // Return value
//   return observation_step;
// }

// INFO: Updated above 2024-08-23
// array[,] vector assemble_observation_step (
//   array[] vector fishing_step,
//   // array[] vector fishing_weight,
//   array[] vector selectivity,
//   vector reporting_step,
//   int K,
//   int L
// ) {
//   // Get dimensions
//   int Years = dims(fishing_step)[1];
//   int S = dims(fishing_step)[2];
//   int T = Years * K;
//   // Declare values
//   array[T, L] vector[S] observation_step;
//   int t = 1;
//   // Populate observation step
//   for (year in 1:Years) {
//     for (k in 1:K) {
//       for (l in 1:L) {
//         observation_step[t, l] = reporting_step
//         .* (1.0 - exp(-fishing_step[year] .* selectivity[l])); // .* fishing_weight[k] * .selectivity[l]
//       }
//       t += 1;
//     }
//   }
//   // Return observation step
//   return observation_step;
// }

array[] vector assemble_fishing_weight (
  array[] vector fishing_weight_transpose
) {
  // Get dimensions
  int W = dims(fishing_weight_transpose)[2];
  int S = dims(fishing_weight_transpose)[1];
  // Declare fishing weight
  array[W] vector[S] fishing_weight;
  // Populate fishing weight
  for (w in 1:W) {
    for (s0 in 1:S) {
      fishing_weight[w, s0] = fishing_weight_transpose[s0, w];
    }
  }
  // Return fishing weight
  return fishing_weight;
}

array[] matrix assemble_movement_rate (
  array[] matrix movement_step,
  int K
) {
  // Get dimensions
  int L = dims(movement_step)[1];
  int S = dims(movement_step)[2];
  // Declare values
  array[L] matrix[S, S] movement_rate;
  // Populate movement rate
  for (l in 1:L) {
    movement_rate[l] = matrix_power(movement_step[l], K);
  }
  // Return movement rate
  return movement_rate;
}

// INFO: Updated above 2024-08-23
// array[] vector assemble_fishing_rate (
//   array[] vector fishing_step,
//   int K
// ) {
//   // Get dimensions
//   int Years = dims(fishing_step)[1];
//   int S = dims(fishing_step)[2];
//   // Declare values
//   array[Years] vector[S] fishing_rate;
//   // Populate fishing rate
//   for (year in 1:Years) {
//     fishing_rate[year] = fishing_step[year] * K;
//   }
//   // Return fishing rate
//   return fishing_rate;
// }

array[] int index_t_to_r (int start, int end) {
  // Declare values
  array[end] int t_to_r;
  int r = 1;
  // Populate index array
  for (t in start:end) {
    t_to_r[t] = r;
    r += 1;
  }
  // Return array
  return t_to_r;
}

real partial_sum_lpmf (
  array[] int index,
  int start,
  int end,
  int T,
  int D,
  int L,
  int S,
  array[,,,,] int tags_transpose,
  array[,] vector tags_released,
  array[,] matrix transition_step,
  array[,] vector observation_step,
  array[] matrix movement_possible,
  real initial_loss_step,
  real tolerance_expected,
  real dispersion
) {
  // Declare index limits
  int R = (end - start + 1); // R stands in for T - 1
  int C = R * D * L * S * S;
  // Declare index arrays
  array[end] int t_to_r = index_t_to_r(start, end);
  // Declare enumeration values
  array[R, D, L] matrix[S, S] abundance;
  array[R, D, L] matrix[S, S] predicted;
  array[C] int observed;
  array[C] real expected;
  // Initialize count
  int count = 0;
  // Populate released abundance
  for (t in start:end) { // Model step
    for (l in 1:L) { // Released size
      abundance[t_to_r[t], 1, l] = diag_matrix(
        tags_released[t, l] * (1 - initial_loss_step)
      );
    }
  }
  // Compute expected recoveries
  for (t in start:end) { // Partial sum index range within released step
    for (d in 2:min(T - t + 1, D)) { // Duration at large
      for (l in 1:L) { // Released size
        // Propagate abundance
        abundance[t_to_r[t], d, l] = abundance[t_to_r[t], d - 1, l]
        * transition_step[t + d - 2, l]; // Previous step
        // Compute predicted
        predicted[t_to_r[t], d, l] = diag_post_multiply(
          abundance[t_to_r[t], d, l],
          observation_step[t + d - 1, l] // Current step
        );
        // Compute vectors
        for (s in 1:S) { // Current region
          for (s0 in 1:S) { // Released region
            if (tags_released[t, l, s0] > 0) { // Were any tags released?
              if (movement_possible[d][s0, s] > 0) {
                // Increment observation count
                count += 1;
                // Populate observed and expected values
                observed[count] = tags_transpose[t, d, l, s, s0]; // Integer
                expected[count] = predicted[t_to_r[t], d, l, s0, s]
                + tolerance_expected; // Real
              } // End if
            } // End if
          } // End s0
        } // End s
      } // End l
    } // End d
  } // End t
  // Return likelihood contribution
  return neg_binomial_2_lupmf(observed[1:count] | expected[1:count],dispersion);
}
