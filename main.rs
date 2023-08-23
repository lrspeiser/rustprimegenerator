#[allow(dead_code)]
fn gen_primes_upto(n: usize) -> Vec<usize> {
    if n < 2 {
        return vec![];
    }

    let mut table = vec![true; n / 2];
    let sqrtn = (n as f64).sqrt().ceil() as usize;

    for i in 2..sqrtn {
        if table[i] {
            for j in (i * i..n).step_by(i) {
                if j % 2 == 0 { continue; }  // Skip even numbers
                table[j / 2] = false;       // Adjusted for half-sized table
            }
        }
    }

    let mut primes: Vec<usize> = vec![2];
    for (i, &is_prime) in table.iter().enumerate().skip(1) {
        if is_prime {
            primes.push(i * 2 + 1);  // Convert index back to odd number
        }
    }
    primes
}

#[allow(dead_code)]
fn gen_primes_upto_segmented(n: usize) -> Vec<usize> {
    if n < 2 {
        return vec![];
    } else if n < 11 {
        return vec![2, 3, 5, 7].into_iter().filter(|&p| p < n).collect();
    }

    let segsize = (n as f64).sqrt().ceil() as usize;
    let baseprimes = gen_primes_upto(segsize);
    let mut primes = baseprimes.clone();

    for segstart in (segsize..n).step_by(segsize) {
        let mut seg = vec![true; segsize];
        let seg_len = seg.len();

        for &bp in &baseprimes {
            let first_multiple = if segstart % bp == 0 {
                segstart
            } else {
                segstart + bp - segstart % bp
            };
            for q in (first_multiple..segstart + segsize).step_by(bp) {
                if q % 2 == 0 { continue; }  // Skip even numbers
                seg[q % seg_len] = false;
            }
        }

        let start = if segstart % 2 == 0 { 1 } else { 0 };
        for i in (start..seg.len()).step_by(2) {
            if seg[i] && (segstart + i) < n {
                primes.push(segstart + i);
            }
        }
    }
    primes
}


fn gen_primes() -> impl Iterator<Item = usize> {
    let mut d = std::collections::HashMap::new();
    let mut q = 2;

    std::iter::from_fn(move || {
        loop {
            if let Some(primes) = d.remove(&q) {
                for &p in &primes {
                    let mut next_composite = p + q;
                    while d.contains_key(&next_composite) {
                        next_composite += p;
                    }
                    d.entry(next_composite).or_insert_with(Vec::new).push(p);
                }
            } else {
                d.insert(q * q, vec![q]);
                let curr_q = q;
                q += 1;
                return Some(curr_q);
            }
            q += 1;
        }
    })
}

// Wheel Factorization
#[allow(dead_code)]
fn gen_primes_upto_wheel(n: usize) -> Vec<usize> {
    if n < 2 {
        return vec![];
    }

    let mut table = vec![true; n + 1];
    let mut primes = vec![2, 3, 5];
    let wheel = [4, 2, 4, 2, 4, 6, 2, 6];

    for &p in &primes {
        for j in (p * p..=n).step_by(p) {
            table[j] = false;
        }
    }

    let mut w = 0;
    let mut num = 7;
    while num <= n {
        if table[num] {
            primes.push(num);
            for j in (num * num..=n).step_by(num) {
                table[j] = false;
            }
        }
        num += wheel[w];
        w = (w + 1) % 8;
    }

    primes
}

// Trial Division
#[allow(dead_code)]
fn gen_primes_upto_trial_division(n: usize) -> Vec<usize> {
    if n < 2 {
        return vec![];
    }

    let mut primes = vec![2];
    for num in 3..=n {
        let sqrtn = (num as f64).sqrt() as usize;
        if primes.iter().take_while(|&&p| p <= sqrtn).all(|p| num % p != 0) {
            primes.push(num);
        }
    }

    primes
}

#[allow(dead_code)]
fn gen_primes_upto_sundaram(n: usize) -> Vec<usize> {
    if n < 2 {
        return vec![];
    }

    let n_half = (n - 1) / 2;
    let mut sieve = vec![true; n_half + 1];

    let limit = ((n_half as f64).sqrt() / 2.0).floor() as usize;
    for i in 1..=limit {
        let mut j = i;
        while i + j + 2 * i * j <= n_half {
            sieve[i + j + 2 * i * j] = false;
            j += 1;
        }
    }

    let mut primes = vec![2];
    for i in 1..=n_half {
        if sieve[i] {
            primes.push(2 * i + 1);
        }
    }

    primes
}


fn main() {
    let n = 100;  // For example, to generate primes up to 1000

    // Measure time for gen_primes
    println!("Using a lazy generation approach with a HashMap:");
    println!("This approach generates prime numbers indefinitely in a lazy manner, i.e., on-demand. It utilizes a HashMap to keep track of the composites and their prime factors, ensuring only the necessary calculations are done. This method is suitable when you need primes on-the-fly without a known upper limit or want to generate primes in real-time scenarios without precomputing a list.");
    let start_gen = std::time::Instant::now();
    let gen = gen_primes();
    let primes_gen: Vec<_> = gen.take_while(|&prime| prime <= n).collect();
    let duration_gen = start_gen.elapsed();
    println!("Primes up to {} using gen_primes: {:?}", n, primes_gen);
    println!("Time taken using gen_primes: {:?}", duration_gen);
    println!("\n");  // Carriage return

    // Measure time for gen_primes_upto_trial_division
    println!("Using the Trial Division approach:");
    println!("Trial Division is a straightforward method that checks the divisibility of a number by primes less than its square root. While simple, it's generally slower than the Sieve of Eratosthenes, especially for larger ranges. However, it's useful when determining the primality of individual numbers without the need for a precomputed list of primes.");
    let start_trial = std::time::Instant::now();
    let primes_trial = gen_primes_upto_trial_division(n);
    let duration_trial = start_trial.elapsed();
    println!("Primes up to {} using gen_primes_upto_trial_division: {:?}", n, primes_trial);
    println!("Time taken using gen_primes_upto_trial_division: {:?}", duration_trial);
    println!("\n");  // Carriage return

    // Measure time for gen_primes_upto_segmented
    println!("Using the Segmented Sieve approach:");
    println!("The Segmented Sieve is a space-optimized version of the Sieve of Eratosthenes. It divides the number range into smaller segments and computes primes in each segment separately, which allows it to generate primes in a range without using memory proportional to the size of the range. This is particularly useful when generating primes in a large range where memory usage is a concern.");
    let start_segmented = std::time::Instant::now();
    let primes_segmented = gen_primes_upto_segmented(n);
    let duration_segmented = start_segmented.elapsed();
    println!("Primes up to {} using gen_primes_upto_segmented: {:?}", n, primes_segmented);
    println!("Time taken using gen_primes_upto_segmented: {:?}", duration_segmented);
    println!("\n");  // Carriage return
  
    // Measure time for gen_primes_upto_wheel
    println!("Using the Wheel Factorization approach:");
    println!("Wheel Factorization optimizes the Sieve of Eratosthenes by skipping over the multiples of the first few primes (e.g., 2, 3, 5). This reduces the number of operations, especially for large numbers, but its effectiveness diminishes as the range grows.");
    let start_wheel = std::time::Instant::now();
    let primes_wheel = gen_primes_upto_wheel(n);
    let duration_wheel = start_wheel.elapsed();
    println!("Primes up to {} using gen_primes_upto_wheel: {:?}", n, primes_wheel);
    println!("Time taken using gen_primes_upto_wheel: {:?}", duration_wheel);
    println!("\n");  // Carriage return

    // Measure time for gen_primes_upto
    println!("Using the basic Sieve of Eratosthenes approach:");
    println!("The Sieve of Eratosthenes is one of the most efficient ways to find all primes smaller than a given number, up to 10 million or so. It works by iteratively marking the multiples of each prime number starting from 2. This method is straightforward and effective for generating a list of primes up to a specified limit.");
    let start_simple = std::time::Instant::now();
    let primes_simple = gen_primes_upto(n);
    let duration_simple = start_simple.elapsed();
    println!("Primes up to {} using gen_primes_upto: {:?}", n, primes_simple);
    println!("Time taken using gen_primes_upto: {:?}", duration_simple);
    println!("\n");  // Carriage return

    // Measure time for gen_primes_upto_sundaram
    println!("Using the Sieve of Sundaram approach:");
    println!("The Sieve of Sundaram is a variant of the Sieve of Eratosthenes. It works by eliminating numbers that can be written in a certain form, reducing the set of numbers to check. The remaining numbers are then transformed to produce a list of primes. This method is simpler in construction than the Sieve of Eratosthenes but can be less efficient for larger ranges.");
    let start_sundaram = std::time::Instant::now();
    let primes_sundaram = gen_primes_upto_sundaram(n);
    let duration_sundaram = start_sundaram.elapsed();
    println!("Primes up to {} using gen_primes_upto_sundaram: {:?}", n, primes_sundaram);
    println!("Time taken using gen_primes_upto_sundaram: {:?}", duration_sundaram);
    println!("\n");  // Carriage return
}
