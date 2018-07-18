#ifndef __MERGETRIMREADS_HPP
#define __MERGETRIMREADS_HPP

using namespace __shedskin__;
namespace __MergeTrimReads__ {

extern str *const_0, *const_1, *const_10, *const_11, *const_12, *const_13, *const_14, *const_15, *const_16, *const_17, *const_18, *const_19, *const_2, *const_20, *const_21, *const_22, *const_23, *const_24, *const_25, *const_26, *const_27, *const_28, *const_29, *const_3, *const_30, *const_31, *const_32, *const_33, *const_34, *const_35, *const_36, *const_37, *const_38, *const_4, *const_5, *const_6, *const_7, *const_8, *const_9;


typedef __ss_int (*lambda0)(str *);
typedef str *(*lambda1)(__ss_int);

extern __ss_bool __3, __4, handle_key, options_allowMissing, options_mergeoverlap, options_onlyoverlap;
extern __ss_int len_key1, len_key2, maxadapter_comp, min_length, offset, options_min_overlap_seqs, options_trimCutoff;
extern tuple2<str *, str *> *keys;
extern str *__name__, *options_adapter_chimera, *table;
extern list<str *> *adapter_chimeras, *options_adapter_F, *options_adapter_S;
extern double cutoff_merge_seqs, cutoff_merge_seqs_early, cutoff_merge_trim, max_prob_N, max_prob_same_channel;


extern __ss_bool  default_1;
extern __ss_bool  default_2;
extern __ss_bool  default_3;
extern list<double> * default_0;
extern str * default_4;
extern str * default_5;
extern str * default_6;

__ss_int edits(str *seq1, str *seq2);
double quality_ident(str *seq1, list<double> *qual1, str *seq2, list<double> *qual2, __ss_int maxcomp);
list<__ss_int> *convert_quality_logprob(str *qualstring);
str *convert_logprob_quality(list<__ss_int> *probs);
list<double> *convert_logprob_prob(list<__ss_int> *lprobs);
str *revcompl(str *seq);
tuple2<str *, __ss_int> *cons_base_prob(str *base1, str *base2, double prob1, double prob2);
tuple2<str *, str *> *check_merge(str *read1, list<double> *qual1, list<__ss_int> *pqual1, str *read2, list<double> *qual2, list<__ss_int> *pqual2);
void *set_options(__ss_int trimcutoff, __ss_bool allowMissing, __ss_bool mergeoverlap, __ss_bool onlyoverlap, __ss_int min_overlap_seqs);
void *set_adapter_sequences(str *forward, str *reverse, str *chimera, __ss_int max_comp);
__ss_bool set_keys(str *key_text);
tuple2<str *, str *> *process_PE(str *read1, str *qual1, str *read2, str *qual2);
tuple2<str *, str *> *overlap_reads(str *read1, str *qual1, str *read2, str *qual2);
tuple2<str *, str *> *consensus_reads(str *read1, str *qual1, str *read2, str *qual2);
tuple2<str *, str *> *process_SR(str *read1, str *qual1);

extern "C" {
PyMODINIT_FUNC initMergeTrimReads(void);

PyMODINIT_FUNC addMergeTrimReads(void);

}
} // module namespace
namespace __shedskin__ { /* XXX */

}
#endif
