use std::cmp::{max, Ordering::Equal};
use std::collections::HashMap;
use std::f64::consts::E;
use mathru::special::gamma::gamma_u;

// LENGTH OF DIRECTIONAL RUNS
pub fn length_of_directional_runs(sequence: &Vec<i8>) -> usize {
    // Step 1
    // Construct the sequence where:
    // if si > si+1 push -1 to s_bar
    // if si <= si+1 push 1 to s_bar
    let s_bar = sequence.windows(2).map(|el| if el[0] <= el[1] { 1 } else { -1 }).collect::<Vec<i8>>();

    // Test statistic T is the length of the longest run in s_bar
    // Compute runs
    let mut max_run = 1;
    let mut curr_run = 1;

    s_bar.windows(2).for_each(|el| {
        // if el[0] != el[1] start the new run
        // else curr_run + 1
        curr_run = if el[0] != el[1] { 1 } else { curr_run + 1 };
        max_run = max(max_run, curr_run);
    });

    println!("max length of run: {}", max_run);
    max_run
}

// INDEPENDENCE FOR NON BINARY DATA
#[derive(PartialEq, Eq, Hash, Debug, Clone, Copy, PartialOrd, Ord)]
pub struct Pair(usize, usize);

#[derive(PartialEq, Debug, Clone, PartialOrd)]
pub struct Bin {
    bin: usize,
    pairs: Vec<Pair>,
    expectation: f64,
    observed_frequency: usize
}

pub fn independence_for_non_binary_data(sequence: &Vec<usize>) -> (f64, usize) {
    // Step 1
    // L = sequence.len()
    // Find the proportion pi of each xi in S
    // number of xi in S / L
    let len = sequence.len();
    let mut props: HashMap<usize, f64> = HashMap::new();

    for el in sequence {
        *props.entry(*el).or_insert(0.0) += 1.0;
    }

    props.iter_mut().for_each(|(_, prop)| *prop /= len as f64);
        
    // Calculate the expected number of occurrences of each possible pair
    // (zi, zj) in S
    // ei,j = pi * pj * L / 2
    let mut pairs: HashMap<Pair, f64> = HashMap::new();

    for el in sequence.windows(2) {
        pairs.entry(Pair(el[0], el[1])).or_insert(props[&el[0]] * props[&el[1]] * (len / 2) as f64);
    }

    // Step 2
    // Allocate the possible (zi, zj pairs, starting from the smallest ei,j
    // into bins such that the expected value of each bin is at least five.
    // The expected value of a bin is equal to the sum of
    // the ei,j values of the pairs that are included in the bin.

    // Convert vector sort by proportion of Pair
    let mut sorted_pairs = Vec::from_iter(pairs);
    sorted_pairs.sort_by(|&(_, ex0), (_, ex1)| ex0.partial_cmp(ex1).unwrap_or(Equal));

    // Bubble sort for organizing repeated proportions
    let mut bucket;
    for i in 0..sorted_pairs.len() - 1 {
        if sorted_pairs[i].1 == sorted_pairs[i+1].1 && sorted_pairs[i].0.0 > sorted_pairs[i+1].0.0 {
            bucket = sorted_pairs[i];
            sorted_pairs[i] =  sorted_pairs[i+1];
            sorted_pairs[i+1] = bucket;
        }
    }

    // After allocating all pairs, if the expected
    // value of the last bin is less than five, merge the last two bins.
    // Let nbin be the number of bins constructed using this procedure.
    let mut bin = 1;

    let mut bins = vec![Bin {
        bin,
        pairs: vec![sorted_pairs[0].0],
        expectation: sorted_pairs[0].1,
        observed_frequency: 0
    }];

    for (pair, expect) in sorted_pairs.iter().skip(1) {
        if let Some(last_expect) = bins.last_mut() {
            if last_expect.expectation >= 5.0 {
                bin+=1;
                bins.push(Bin {bin, pairs: vec![*pair], expectation: *expect, observed_frequency: 0 });
            }
            else {
                last_expect.pairs.push(*pair);
                last_expect.expectation += *expect;
            }
        }
    }

    // Step 3 Chi-square test
    // Let o be a list of nbin counts, each initialized to 0. For j=1 to L-1:
    // a. If the pair (sj, sj+1) is in bin i, increment oi by 1.
    // b. Let j = j+2.
    for el in sequence.windows(2).step_by(2) {
        for bin in &mut bins {
            if bin.pairs.contains(&Pair(el[0], el[1])) {
                    bin.observed_frequency +=1;
            }
        }
    }

    // Compute test statistic
    // i = 1 to number of bins
    // (oi - expectation)^2 / expectation
    let mut t:f64 = bins.iter().fold(0.0, |acc, bin| acc + ((bin.observed_frequency as f64 - bin.expectation).powi(2) / bin.expectation));
    t = (t * 100.0).round() / 100.0;

    // degrees of freedom = nbin - k
    let dof = bins.len() - props.len();

    println!("t: {}", t);
    println!("degrees of freedom: {}", dof);
        
    (t, dof)
    }

// THE COLLISION ESTIMATE
pub fn collision_estimate(bits: &Vec<i8>) -> f64 {

    // Step 1
    // index = 1
    // t = vector of t's
    let mut index = 1;
    let mut t: Vec<usize> = Vec::new();

    // Step 2
    // Beginning with S index, step through the input until any observed value is repeated
    // find the smallest j such that si = sj, for some i with index ≤ i < j
    // S = bits
    // Step 3 + 4
    // Repeat steps 2-3 until the end of the dataset is reached
    while index < bits.len() {
        // 11 00
        if bits[index -1] == bits[index] {
            t.push(2);
            index += 2;
        }
        // 101 011 
        else if index < bits.len() -1 {
            t.push(3); 
            index += 3;
        }
        else {
            break;
        }
    }

    // Step 5
    // Calculate the sample mean X, and the sample standard deviation
    // t_len = length of vector t
    let t_len: usize = t.len();
    // t_sum = sum of t's elements
    let t_sum: usize = t.iter().sum();
    // x_bar = sample mean
    let mean = t_sum as f64 / t_len as f64;
    // sq_sum = temp sum for standard deviation
    let sq_sum: f64 = t.iter().map(|el| (*el as f64 - mean).powi(2)).sum();


    // Step 6
    // Compute the lower-bound of the confidence interval for the mean,
    // based on a normal distribution with a confidence level of 99 %

    // sigma_hat = standard deviation
    let sigma_hat = (sq_sum / (t_len - 1) as f64).sqrt();
    // x_bar_prime = lower-bound of the confidence interval for the mean
    let x_bar_prime = mean - 2.576 * (sigma_hat / (t_len as f64).sqrt());

    // Step 7
    // Solve the p by binary search between 0.5 & 1
    let mut p = 0.0;
    let mut min_p = 0.5;
    let mut mid_p = 0.75;
    let mut max_p = 1.0;
    let mut possible_x_bar_prime;

        for _ in 0..10_000 {
            // get possible_x_bar_prime with mid_p
            possible_x_bar_prime = f_x_bar_prime(mid_p);

            // p is greater than we have 
            if possible_x_bar_prime > x_bar_prime {
                min_p = mid_p;
            }

            // p is smaller than we have
            else if possible_x_bar_prime < x_bar_prime {
                max_p = mid_p;
            }

            // correct p found
            else if possible_x_bar_prime == x_bar_prime {
                p = mid_p;
                break;
            }
            // update the mid_p for next iteration
            mid_p = (min_p + max_p) / 2.0;
        }

        // Step 8
        // If the binary search yields a solution
        // min_entropy = -log2(p)
        // If the search does not yield a solution
        // min_entropy = log2(2) = 1
        let min_entropy =  if t_len != 0 && p < 1.0 && p > 0.5 { -p.log2() } else { 1.0 };
        
        println!("t: {:?}", t);
        println!("length of t: {}",t_len);
        println!("mean: {:.4}", mean );
        println!("sigma hat: {:.4}", sigma_hat);
        println!("x bar prime: {:.4}",  x_bar_prime);
        println!("p: {:.4}", p);
        println!("min entropy: {:.4}", min_entropy);

        min_entropy

}

fn f_x_bar_prime(p: f64) -> f64 {
    let q = 1.0 - p;
    let z = 1.0 / q;

    let incomplete_gamma_result = gamma_u(3.0, z) * (z.powi(-3)) * (E.powf(z));

    (p * (q.powf(-2.0))) * (1.0 + 0.5 * ((1.0 / p) - (1.0 / q))) * incomplete_gamma_result -
    (p * (q.powf(-1.0))) * 0.5 * ((1.0 / p) - (1.0 / q))
}

#[cfg(test)]
pub mod test {
    use super::*;

    #[test]
    pub fn test_length_of_directional_runs() {
        let input: Vec<i8> = vec![2, 2, 2, 5, 7, 7, 9, 3, 1, 4, 4];
        let longest_run = length_of_directional_runs(&input);

        assert_eq!(longest_run, 6);
    }

    #[test]
    pub fn test_independence_for_non_binary_data() {
        let sequence: Vec<usize> = vec![2, 2, 3, 1, 3, 2, 3, 2, 1, 3, 1, 1, 2, 3, 1, 1, 2, 2, 2, 3, 3, 2, 3, 2, 3, 1, 2, 2, 3, 3,
        2, 2, 2, 1, 3, 3, 3, 2, 3, 2, 1, 3, 2, 3, 1, 2, 2, 3, 1, 1, 3, 2, 3, 2, 3, 1, 2, 2, 3, 3, 2, 2, 2, 1, 3, 3, 3, 2, 3,
        2, 1, 2, 2, 3, 3, 3, 2, 3, 2, 1, 2, 2, 2, 1, 3, 3, 3, 2, 3, 2, 1, 3, 2, 3, 1, 2, 2, 3, 1, 1];
        let (t, d) = independence_for_non_binary_data(&sequence);

        assert_eq!(t, 3.46);
        assert_eq!(d, 3);
    }

    #[test]
    pub fn test_collision_estimate() {
        let bits: Vec<i8> = vec![1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1,
        0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0];
        let min_entropy = collision_estimate(&bits);

        assert_eq!(min_entropy, 0.448339587791084);
    }
}