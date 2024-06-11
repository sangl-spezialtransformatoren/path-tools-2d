use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::hash::{Hash, Hasher};
use std::ops::Bound::Included;

use crate::geom::Point;
use crate::geom::Segment;
use crate::geom::SegmentTrait;
use crate::helpers::ApproximateEq;

#[derive(Clone)]
struct IndexedSegment {
    index: usize,
    segment: Segment,
}

impl IndexedSegment {
    fn new(index: usize, segment: Segment) -> Self {
        IndexedSegment { index, segment }
    }

    fn reversed(&self) -> Self {
        IndexedSegment {
            index: self.index,
            segment: self.segment.reversed(),
        }
    }
}

fn ord_float(a: f64, b: f64) -> Ordering {
    a.partial_cmp(&b).expect("Invalid value")
}

fn ord_points(a: &Point, b: &Point) -> Ordering {
    if a.x.approximate_eq(b.x) {
        ord_float(a.y, b.y)
    } else {
        ord_float(a.x, b.x)
    }
}

impl Eq for IndexedSegment {}

impl PartialEq<Self> for IndexedSegment {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index
    }
}

impl PartialOrd<Self> for IndexedSegment {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for IndexedSegment {
    fn cmp(&self, other: &Self) -> Ordering {
        match ord_points(&self.segment.start(), &other.segment.start()) {
            Ordering::Equal => self.index.cmp(&other.index),
            order => order,
        }
    }
}

fn find_and_remove(
    primary_tree: &mut BTreeSet<IndexedSegment>,
    secondary_tree: &mut BTreeSet<IndexedSegment>,
    segment: &Segment,
) -> Option<Segment> {
    let next_element = primary_tree
        .range((
            Included(IndexedSegment::new(usize::MIN, segment.clone())),
            Included(IndexedSegment::new(usize::MAX, segment.clone())),
        ))
        .next()
        .cloned();

    if let Some(item) = next_element {
        primary_tree.remove(&item);
        secondary_tree.remove(&item.reversed());
        Some(item.segment)
    } else {
        None
    }
}

pub fn combine_segments(segments: &[Segment]) -> Vec<Vec<Segment>> {
    let mut forward_tree: BTreeSet<IndexedSegment> = BTreeSet::new();
    let mut backward_tree: BTreeSet<IndexedSegment> = BTreeSet::new();

    for (i, segment) in segments.iter().enumerate() {
        let segment = segment.clone();
        forward_tree.insert(IndexedSegment::new(i, segment));
        backward_tree.insert(IndexedSegment::new(i, segment.reversed()));
    }

    let mut result: Vec<Vec<Segment>> = vec![];

    while let Some(IndexedSegment { index: i, segment }) = forward_tree.pop_first() {
        backward_tree.remove(&IndexedSegment::new(i, segment.reversed()));
        let mut current_strain = vec![segment];

        loop {
            let current_segment = current_strain
                .last()
                .expect("Current strain must have elements");

            if let Some(next_segment) = find_and_remove(
                &mut forward_tree,
                &mut backward_tree,
                &current_segment.reversed(),
            ) {
                current_strain.push(next_segment);
                continue;
            }

            let current_segment = current_strain
                .first()
                .expect("Current strain must have elements");

            if let Some(next_segment) =
                find_and_remove(&mut backward_tree, &mut forward_tree, current_segment)
            {
                current_strain.insert(0, next_segment.reversed());
                continue;
            }

            break;
        }

        result.push(current_strain);
    }

    result
}
