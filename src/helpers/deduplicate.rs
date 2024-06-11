pub fn deduplicate<T>(values: &[T], is_equal: fn(&T, &T) -> bool) -> Vec<&T> {
    let mut deduplicated_values = Vec::with_capacity(values.len());

    'outer: for (i, value) in values.iter().enumerate() {
        for (j, other_value) in values[0..i].iter().enumerate() {
            if i != j && is_equal(value, other_value) {
                continue 'outer;
            }
        }
        deduplicated_values.push(value);
    }
    deduplicated_values
}
