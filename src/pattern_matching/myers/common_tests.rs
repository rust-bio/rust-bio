macro_rules! impl_tests {
    ($mod_:ident, $bitvec:ty, $dist_type:ty, $builder_method:ident) => {
        use crate::alignment::AlignmentOperation::*;
        use crate::alignment::{Alignment, AlignmentMode};
        use crate::pattern_matching::myers::MyersBuilder;
        use itertools::Itertools;
        use $mod_::Myers;

        #[test]
        fn test_find_all_end() {
            let text = "ACCGTGGATGAGCGCCATAG".to_string();
            let patt = "------GATGAGCGT-----".replace('-', "");
            let myers = Myers::<$bitvec>::new(patt.as_bytes());
            let occ = myers.find_all_end(text.as_bytes(), 1).collect_vec();
            assert_eq!(occ, [(13, 1), (14, 1)]);
        }

        #[test]
        fn test_distance() {
            let text = b"TGAGCNTA";
            let patt = b"TGAGCGT";

            let myers = Myers::<$bitvec>::new(patt);
            assert_eq!(myers.distance(text), 1);

            let myers_wildcard = MyersBuilder::new().text_wildcard(b'N').build_64(patt);
            assert_eq!(myers_wildcard.distance(text), 0);
        }

        #[test]
        fn test_full_position() {
            let text = "CAGACATCTT".to_string();
            let patt = "-AGA------".replace('-', "");

            let mut myers = Myers::<$bitvec>::new(patt.as_bytes());
            let matches: Vec<_> = myers.find_all(text.as_bytes(), 1).collect();
            assert_eq!(&matches, &[(1, 3, 1), (1, 4, 0), (1, 5, 1), (3, 6, 1)]);
        }

        #[test]
        fn test_traceback_path() {
            let text = "TCAGACAT-CTT".replace('-', "");
            let patt = "TC-GACGTGCT".replace('-', "");

            let mut myers = Myers::<$bitvec>::new(patt.as_bytes());
            let mut matches = myers.find_all(text.as_bytes(), 3);
            let mut aln = vec![];
            assert_eq!(matches.next_path(&mut aln).unwrap(), (0, 10, 3));
            assert_eq!(
                aln,
                &[Match, Match, Del, Match, Match, Match, Subst, Match, Ins, Match, Match]
            );
        }

        #[test]
        fn test_traceback_path2() {
            let text = "TCAG--CAGATGGAGCTC".replace('-', "");
            let patt = "TCAGAGCAG---------".replace('-', "");

            let mut myers = Myers::<$bitvec>::new(patt.as_bytes());
            let mut matches = myers.find_all(text.as_bytes(), 2);
            let mut aln = vec![];
            assert_eq!(matches.next_path(&mut aln).unwrap(), (0, 7, 2));
            assert_eq!(
                aln,
                &[Match, Match, Match, Match, Ins, Ins, Match, Match, Match]
            );
        }

        #[test]
        fn test_alignment() {
            let text = "GGTCCTGAGGGATTA".to_string();
            let patt = "--TCCT-AGGGA---".replace('-', "");

            let mut myers = Myers::<$bitvec>::new(patt.as_bytes());
            let expected = Alignment {
                score: 1,
                xstart: 0,
                xend: 9,
                xlen: 9,
                ystart: 2,
                yend: 12,
                ylen: 15,
                operations: vec![
                    Match, Match, Match, Match, Del, Match, Match, Match, Match, Match,
                ],
                mode: AlignmentMode::Semiglobal,
            };

            let mut aln = Alignment::default();
            {
                let mut matches = myers.find_all(text.as_bytes(), 1);
                assert!(matches.next_alignment(&mut aln));
                assert_eq!(&aln, &expected);

                aln.score = -1;
                matches.alignment(&mut aln);
                assert_eq!(&aln, &expected);
            }
            // Lazy API
            aln.score = -1;
            let end = expected.yend - 1;
            let mut lazy_matches = myers.find_all_lazy(text.as_bytes(), 1);
            assert!(!lazy_matches.alignment_at(end, &mut aln));
            // now search positions
            aln.score = -1;
            assert_eq!(
                lazy_matches.next(),
                Some((end, expected.score as $dist_type))
            );
            assert!(lazy_matches.alignment_at(end, &mut aln));
            assert_eq!(&aln, &expected);
        }

        #[test]
        fn test_position_cmp() {
            // same as position_at, but 0-based positions from
            let text = "CAGACATCTT".to_string();
            let patt = "---AGA----".replace('-', "");
            let text = text.as_bytes();

            let starts_exp = [1, 1, 1, 3];
            let end_dist_exp = [(2, 1), (3, 0), (4, 1), (5, 1)];

            let mut myers = Myers::<$bitvec>::new(patt.as_bytes());

            // standard iterator with 0-based ends
            let end_dist: Vec<_> = myers.find_all_end(text, 1).collect();
            assert_eq!(&end_dist, &end_dist_exp);

            // iterator over full ranges where ends are + 1
            let full_hits: Vec<_> = myers.find_all(text, 1).collect();

            // allows to retrive starting position later
            let mut lazy_matches = myers.find_all_lazy(text, 1);

            // compare with each other and lazy matches
            for ((&start, (end, dist)), (f_start, f_end, f_dist)) in
                starts_exp.iter().zip(end_dist).zip(full_hits)
            {
                assert_eq!(start, f_start);
                assert_eq!(dist, f_dist);
                assert_eq!(end + 1, f_end);

                // lazy API
                let (lazy_end, lazy_dist) = lazy_matches.next().unwrap();
                assert_eq!(end, lazy_end);
                assert_eq!(dist, lazy_dist);
                assert_eq!(lazy_matches.hit_at(end), Some((start, dist)));
                // For positions above, information is not (yet) available
                assert_eq!(lazy_matches.hit_at(end + 1), None);
            }
        }

        #[test]
        fn test_path_at() {
            let text = "CAGACATCTT".to_string();
            let patt = "---AGA----".replace('-', "");

            let mut myers = Myers::<$bitvec>::new(patt.as_bytes());
            let mut matches = myers.find_all_lazy(text.as_bytes(), 1);

            let expected = &[Match, Match, Ins];
            let mut path = vec![];

            // search first hit
            assert_eq!(matches.next(), Some((2, 1)));

            // retrieve first hit at 0-based end position (2)
            assert_eq!(matches.hit_at(2), Some((1, 1)));
            assert_eq!(matches.path_at(2, &mut path), Some((1, 1)));
            assert_eq!(&path, expected);

            // hit out of range
            path.clear();
            assert!(matches.path_at(3, &mut path).is_none());
            assert!(path.is_empty());

            // now search the next hit
            assert_eq!(matches.next(), Some((3, 0)));
            // position 3 is now searched -> path can be retrieved
            assert_eq!(matches.path_at(3, &mut path), Some((1, 0)));
            assert_eq!(&path, &[Match, Match, Match])
        }

        #[test]
        fn test_shorter() {
            let text = "-ATG".replace('-', "");
            let patt = "CATGC".to_string();

            let mut myers = Myers::<$bitvec>::new(patt.as_bytes());
            let mut matches = myers.find_all(text.as_bytes(), 2);
            let mut aln = vec![];
            assert_eq!(matches.next_path(&mut aln).unwrap(), (0, 3, 2));
            assert_eq!(aln, &[Ins, Match, Match, Match, Ins]);
        }

        #[test]
        fn test_long_shorter() {
            let text =
                "C--------CACGCGTGGGTCCTGAGGGAGCTCGTCGGTGTGGGGTTCGGGGGGGTTTGT".replace('-', "");
            let patt = "CGGGGTGTGCACGCGTGGGTCCTGAGGGAGCTCGTCGGTGTGGGGTTCGGGGGGGTTTGT".to_string();

            let mut myers = Myers::<$bitvec>::new(patt.as_bytes());
            let mut matches = myers.find_all(text.as_bytes(), 8);
            assert_eq!(matches.next().unwrap(), (0, 52, 8));
        }

        #[test]
        fn test_ambig() {
            let text = b"TGABCNTR";
            let patt = b"TGRRCGTR";
            //                x  x
            // Matching is asymmetric here (A matches R and G matches N, but the reverse is not true)

            let myers = MyersBuilder::new().ambig(b'R', b"AG").build_64(patt);
            assert_eq!(myers.distance(text), 2);
        }

        #[test]
        fn test_longest_possible() {
            let text = b"CCACGCGT";

            let mut myers: Myers<u8> = Myers::new(text);
            assert_eq!(myers.find_all(text, 0).next(), Some((0, 8, 0)));
        }

        #[test]
        fn test_large_dist() {
            use std::iter::repeat;

            let pattern: Vec<_> = repeat(b'T').take(64).collect();
            let text: Vec<_> = repeat(b'A').take(64).collect();

            let mut myers = Myers::<u64>::new(&pattern);
            let max_dist = myers
                .find_all_end(&text, 64)
                .max_by_key(|&(_, dist)| dist)
                .unwrap()
                .1;
            assert_eq!(max_dist, 64);

            let max_dist = myers
                .find_all(&text, 64)
                .max_by_key(|&(_, _, dist)| dist)
                .unwrap()
                .2;
            assert_eq!(max_dist, 64);
        }

        // macro end
    };
}
