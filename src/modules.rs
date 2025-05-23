use rayon::prelude::*;

// get NX
pub fn get_nx(numbers: &mut [i64], fraction: f32) -> i64 {
    numbers.sort_unstable();

    // half of the bases
    let halfsum = numbers.par_iter().sum::<i64>() as f32 * fraction; // f32 * f32

    // cumsum of the sorted vector
    let cumsum = numbers
        .iter()
        .scan(0, |sum, i| {
            *sum += i;
            Some(*sum)
        })
        .collect::<Vec<_>>();
    let n50_index = cumsum.par_iter().position_first(|&x| x > halfsum as i64).unwrap();

    numbers[n50_index]
}

// get number of bases with q >= value
pub fn get_qual_bases(q: &[u8], qx: u8) -> i64 {
    let mut n = 0;
    for item in q {
        if *item >= qx {
            n += 1
        }
    }
    n
}

pub fn get_gc_bases(seq: &[u8]) -> u64 {
    let mut n: u64 = 0;
    for s in seq {
        if matches!(s, &103u8 | &71u8 |  &99u8 | &67u8) { //G, g, C or c
            n += 1;
        }
    }
    n
}

// to get mean of q scores from a record - first convert to prob, calc mean, then back to phred
// this fn reads phred and converts to probs and returns their sum
//
pub fn qscore_mean(q: &[u8]) -> u8 {
    let mut qprob_sum = 0.0;
    let mut len = 0;
    for &item in q {
        len += 1;
        let phred = *&item as f32 - 33.0;
        let prob = 10.0_f32.powf(-phred / 10.0);
        qprob_sum += prob
    }
    let mean_prob = qprob_sum / len as f32;
    
    (-10.0 * mean_prob.log10()) as u8    
}

pub fn median(numbers: &mut Vec<u8>) -> u8 {
    numbers.sort_unstable();

    let mid = numbers.len() / 2;
    if numbers.len() % 2 == 0 {
        return (numbers[mid - 1] + numbers[mid]) / 2
    } else {
        return numbers[mid]
    }
}

