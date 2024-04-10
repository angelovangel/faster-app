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
