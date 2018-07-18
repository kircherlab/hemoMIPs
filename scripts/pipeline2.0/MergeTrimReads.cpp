#include "builtin.hpp"
#include "os/path.hpp"
#include "os/__init__.hpp"
#include "string.hpp"
#include "sys.hpp"
#include "math.hpp"
#include "stat.hpp"
#include "MergeTrimReads.hpp"

/**

Merge/Adapter trim reads stored in BAM

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *02.01.2012
:Type: module for C cross-compilation

*/

namespace __MergeTrimReads__ {

str *const_0, *const_1, *const_10, *const_11, *const_12, *const_13, *const_14, *const_15, *const_16, *const_17, *const_18, *const_19, *const_2, *const_20, *const_21, *const_22, *const_23, *const_24, *const_25, *const_26, *const_27, *const_28, *const_29, *const_3, *const_30, *const_31, *const_32, *const_33, *const_34, *const_35, *const_36, *const_37, *const_38, *const_4, *const_5, *const_6, *const_7, *const_8, *const_9;


list<str *> *adapter_chimeras, *options_adapter_F, *options_adapter_S;
tuple2<str *, str *> *keys;
double cutoff_merge_seqs, cutoff_merge_seqs_early, cutoff_merge_trim, max_prob_N, max_prob_same_channel;
__ss_bool __3, __4, handle_key, options_allowMissing, options_mergeoverlap, options_onlyoverlap;
str *__name__, *options_adapter_chimera, *table;
__ss_int len_key1, len_key2, maxadapter_comp, min_length, offset, options_min_overlap_seqs, options_trimCutoff;


__ss_bool  default_1;
__ss_bool  default_2;
__ss_bool  default_3;
list<double> * default_0;
str * default_4;
str * default_5;
str * default_6;

static inline __ss_int __lambda0__(str *x);
static inline str *__lambda1__(__ss_int x);

static inline __ss_int __lambda0__(str *x) {
    
    return (ord(x)-offset);
}

static inline str *__lambda1__(__ss_int x) {
    
    return chr(___max(2, 0, 33, ___min(2, 0, 126, (x+offset))));
}

__ss_int edits(str *seq1, str *seq2) {
    __ss_int __5, __6, dist, lmax, lmin, pos;

    lmin = ___min(2, 0, len(seq1), len(seq2));
    lmax = ___max(2, 0, len(seq1), len(seq2));
    dist = (lmax-lmin);

    FAST_FOR(pos,0,lmin,1,5,6)
        if (__ne(seq1->__getfast__(pos), seq2->__getfast__(pos))) {
            dist = (dist+1);
        }
    END_FOR

    return dist;
}

double quality_ident(str *seq1, list<double> *qual1, str *seq2, list<double> *qual2, __ss_int maxcomp) {
    str *hseq1;
    __ss_bool __10, __11, __12, __13, __16, __17, __18, __19, __20, __23, __24, __25, __26, __27, __7, __8, __9;
    __ss_int __14, __15, __21, __22, pos;
    double ident, res, val;
    list<double> *hqual1;

    ident = 0.0;
    res = 0.0;
    if (((len(qual2)!=0) and (len(seq2)<len(seq1)))) {
        hseq1 = seq1;
        hqual1 = qual1;
        qual1 = qual2;
        seq1 = seq2;
        qual2 = hqual1;
        seq2 = hseq1;
    }
    if (((maxcomp<0) or (maxcomp>len(seq1)))) {
        maxcomp = len(seq1);
    }
    if (((len(seq1)>=maxcomp) and (len(seq2)>=maxcomp) and (len(qual1)>=maxcomp))) {
        if ((len(qual2)==0)) {

            FAST_FOR(pos,0,maxcomp,1,14,15)
                if ((__ne(seq1->__getfast__(pos), seq2->__getfast__(pos)) and __ne(seq1->__getfast__(pos), const_0) and __ne(seq2->__getfast__(pos), const_0))) {
                    ident = (ident+(1.0-qual1->__getfast__(pos)));
                }
                else if ((__eq(seq1->__getfast__(pos), const_0) or __eq(seq2->__getfast__(pos), const_0))) {
                    maxcomp = (maxcomp-1);
                }
                else {
                    ident = (ident+1.0);
                }
            END_FOR

        }
        else {

            FAST_FOR(pos,0,maxcomp,1,21,22)
                if ((__ne(seq1->__getfast__(pos), seq2->__getfast__(pos)) and __ne(seq1->__getfast__(pos), const_0) and __ne(seq2->__getfast__(pos), const_0))) {
                    val = (((qual1->__getfast__(pos)<qual2->__getfast__(pos)))?(qual1->__getfast__(pos)):(qual2->__getfast__(pos)));
                    ident = (ident+(1.0-val));
                }
                else if ((__eq(seq1->__getfast__(pos), const_0) or __eq(seq2->__getfast__(pos), const_0))) {
                    maxcomp = (maxcomp-1);
                }
                else {
                    ident = (ident+1.0);
                }
            END_FOR

        }
    }
    else {
        (__sys__::__ss_stderr)->write(const_1);
        print2(NULL,0,5, ___box(___bool((len(seq1)>=maxcomp))), ___box(___bool((len(seq2)>=maxcomp))), ___box(___bool((len(qual1)>=maxcomp))), ___box(len(seq2)), ___box(maxcomp));
        __sys__::__ss_exit(0);
        return 0.0;
    }
    if ((maxcomp>0)) {
        res = (ident/__float(maxcomp));
    }
    return res;
}

list<__ss_int> *convert_quality_logprob(str *qualstring) {
    
    return map(1, __lambda0__, qualstring);
}

str *convert_logprob_quality(list<__ss_int> *probs) {
    
    return (const_2)->join(map(1, __lambda1__, probs));
}

list<double> *convert_logprob_prob(list<__ss_int> *lprobs) {
    __iter<__ss_int> *__29;
    double val;
    list<__ss_int> *__28;
    list<double> *res;
    list<__ss_int>::for_in_loop __31;
    __ss_int __30, prob;

    res = (new list<double>());

    FOR_IN(prob,lprobs,28,30,31)
        val = (1.0-__math__::pow(10.0, (prob/(-10.0))));
        if ((val>max_prob_N)) {
            res->append(val);
        }
        else {
            res->append(max_prob_N);
        }
    END_FOR

    return res;
}

str *revcompl(str *seq) {
    
    return (const_2)->join((seq->translate(table))->__slice__(4, 0, 0, (-1)));
}

tuple2<str *, __ss_int> *cons_base_prob(str *base1, str *base2, double prob1, double prob2) {
    list<str *> *__36, *bases;
    list<str *>::for_in_loop __39;
    double aprob1, aprob2, help, lprob1, lprob2, thelp, total_prob;
    __iter<str *> *__33, *__37, *__41;
    str::for_in_loop __35, __43;
    str *__32, *__40, *call_base, *char_call, *char_elem;
    __ss_int __34, __38, __42, call_quality, hqual, val;

    aprob1 = (__math__::log10((1.0-prob1))-__math__::log10(3.0));
    lprob1 = __math__::log10(prob1);
    lprob2 = __math__::log10(prob2);
    aprob2 = (__math__::log10((1.0-prob2))-__math__::log10(3.0));
    bases = (new list<str *>(2,base1,base2));
    if (__eq(base1, base2)) {
        bases = (new list<str *>(1,base1));
    }
    total_prob = 0.0;

    FOR_IN(char_elem,const_3,32,34,35)
        help = 0.0;
        if (__eq(base1, char_elem)) {
            help = (help+lprob1);
        }
        else {
            help = (help+aprob1);
        }
        if (__eq(base2, char_elem)) {
            help = (help+lprob2);
        }
        else {
            help = (help+aprob2);
        }
        total_prob = (total_prob+__power(10, help));
    END_FOR

    total_prob = __math__::log10(total_prob);
    call_base = const_4;
    call_quality = (-1);

    FOR_IN(char_call,bases,36,38,39)
        thelp = 0.0;

        FOR_IN(char_elem,const_3,40,42,43)
            if (__ne(char_elem, char_call)) {
                help = 0.0;
                if (__eq(base1, char_elem)) {
                    help = (help+lprob1);
                }
                else {
                    help = (help+aprob1);
                }
                if (__eq(base2, char_elem)) {
                    help = (help+lprob2);
                }
                else {
                    help = (help+aprob2);
                }
                thelp = (thelp+__power(10, help));
            }
        END_FOR

        thelp = __math__::log10(thelp);
        val = __int(___round(((-10.0)*(thelp-total_prob))));
        hqual = (((val>60))?(60):(val));
        if ((hqual>call_quality)) {
            call_base = char_call;
            call_quality = hqual;
        }
    END_FOR

    return (new tuple2<str *, __ss_int>(2,call_base,call_quality));
}

tuple2<str *, str *> *check_merge(str *read1, list<double> *qual1, list<__ss_int> *pqual1, str *read2, list<double> *qual2, list<__ss_int> *pqual2) {
    list<str *> *new_lseq;
    tuple2<str *, __ss_int> *__63;
    __ss_bool __44, __45, __50, __51, __52, __53, __54, __55, __56, __57, __58, __59, __60, __61, __62;
    list<__ss_int> *new_qual;
    str *nbase, *new_seq;
    __ss_int __46, __47, __48, __49, lread1, lread2, nqual, pos, spos;
    double oident;

    new_seq = const_2;
    new_qual = (new list<__ss_int>());
    lread1 = len(read1);
    lread2 = len(read2);
    if (((lread1>0) and (lread2>0))) {
        oident = quality_ident(read2, qual2, read1, qual1, (-1));
        if ((oident>cutoff_merge_seqs)) {
            spos = (-1);
            new_lseq = (new list<str *>(read1));
            new_qual = pqual1;

            FAST_FOR(pos,lread1,lread2,1,46,47)
                new_lseq->append(read2->__getfast__(pos));
                new_qual->append(pqual2->__getfast__(pos));
            END_FOR


            FAST_FOR(pos,0,lread2,1,48,49)
                if (((spos==(-1)) and __ne(read2->__getfast__(pos), const_0) and __ne(read1->__getfast__(pos), const_0))) {
                    spos = pos;
                }
                if ((__eq(new_lseq->__getfast__(pos), const_0) or (__eq(new_lseq->__getfast__(pos), const_4) and __ne(read2->__getfast__(pos), const_4) and __ne(read2->__getfast__(pos), const_0)))) {
                    new_lseq->__setitem__(pos, read2->__getfast__(pos));
                    new_qual->__setitem__(pos, pqual2->__getfast__(pos));
                }
                else if (((pos<lread1) and __ne(read1->__getfast__(pos), const_4) and __ne(read2->__getfast__(pos), const_4) and __ne(read1->__getfast__(pos), const_0) and __ne(read2->__getfast__(pos), const_0))) {
                    __63 = cons_base_prob(read1->__getfast__(pos), read2->__getfast__(pos), qual1->__getfast__(pos), qual2->__getfast__(pos));
                    nbase = __63->__getfirst__();
                    nqual = __63->__getsecond__();
                    new_qual->__setitem__(pos, nqual);
                    new_lseq->__setitem__(pos, nbase);
                }
            END_FOR

            if ((spos==(-1))) {
                spos = 0;
            }
            if (options_onlyoverlap) {
                new_seq = (const_2)->join(new_lseq->__slice__(3, spos, (lread1+1), 0));
                new_qual = new_qual->__slice__(3, spos, (lread1+1), 0);
            }
            else {
                new_seq = (const_2)->join(new_lseq);
            }
        }
    }
    return (new tuple2<str *, str *>(2,new_seq,convert_logprob_quality(new_qual)));
}

void *set_options(__ss_int trimcutoff, __ss_bool allowMissing, __ss_bool mergeoverlap, __ss_bool onlyoverlap, __ss_int min_overlap_seqs) {
    
    options_trimCutoff = trimcutoff;
    options_allowMissing = allowMissing;
    options_mergeoverlap = mergeoverlap;
    options_onlyoverlap = onlyoverlap;
    options_min_overlap_seqs = min_overlap_seqs;
    return NULL;
}

void *set_adapter_sequences(str *forward, str *reverse, str *chimera, __ss_int max_comp) {
    list<str *> *__64, *__72, *__80, *__84;
    __iter<__ss_int> *__75;
    list<str *>::for_in_loop __67;
    __ss_bool __70, __71;
    str *adapter_F;
    __iter<str *> *__65;
    list<__ss_int>::for_in_loop __77;
    list<__ss_int> *__74, *helpremove;
    __ss_int __66, __68, __69, __73, __76, __78, __79, __81, __82, __83, __85, i, ind;

    options_adapter_F = (forward->upper())->split(const_5);
    options_adapter_S = (reverse->upper())->split(const_5);
    options_adapter_chimera = chimera->upper();
    adapter_chimeras = options_adapter_chimera->split(const_5);

    FOR_IN(adapter_F,options_adapter_F,64,66,67)
        if ((!(adapter_chimeras)->__contains__(adapter_F))) {
            adapter_chimeras->append(adapter_F);
        }
    END_FOR

    helpremove = (new list<__ss_int>());

    FAST_FOR(ind,0,len(adapter_chimeras),1,68,69)
        adapter_chimeras->__setitem__(ind, (adapter_chimeras->__getfast__(ind))->strip());
        if (((len(adapter_chimeras->__getfast__(ind))>0) and (len(adapter_chimeras->__getfast__(ind))<(maxadapter_comp+1)))) {

            while ((len(adapter_chimeras->__getfast__(ind))<(maxadapter_comp+1))) {
                __72 = adapter_chimeras;
                __73 = ind;
                __72->__setitem__(__73, (__72->__getfast__(__73))->__iadd__(const_0));
            }
        }
        else if ((len(adapter_chimeras->__getfast__(ind))==0)) {
            helpremove->append(ind);
        }
    END_FOR


    FOR_IN(ind,helpremove->__slice__(4, 0, 0, (-1)),74,76,77)
        adapter_chimeras->pop(ind);
    END_FOR


    FAST_FOR(i,0,len(options_adapter_F),1,78,79)
        if ((len(options_adapter_F->__getfast__(i))<maxadapter_comp)) {

            while ((len(options_adapter_F->__getfast__(i))<maxadapter_comp)) {
                __80 = options_adapter_F;
                __81 = i;
                __80->__setitem__(__81, (__80->__getfast__(__81))->__iadd__(const_0));
            }
        }
    END_FOR


    FAST_FOR(i,0,len(options_adapter_S),1,82,83)
        if ((len(options_adapter_S->__getfast__(i))<maxadapter_comp)) {

            while ((len(options_adapter_S->__getfast__(i))<maxadapter_comp)) {
                __84 = options_adapter_S;
                __85 = i;
                __84->__setitem__(__85, (__84->__getfast__(__85))->__iadd__(const_0));
            }
        }
    END_FOR

    return NULL;
}

__ss_bool set_keys(str *key_text) {
    __ss_bool __86, __87;

    key_text = key_text->upper();
    if ((key_text->count(const_5)==1)) {
        keys = (new tuple2<str *, str *>(2,(key_text->split(const_5))->__getfast__(0),(key_text->split(const_5))->__getfast__(1)));
    }
    else if ((key_text->count(const_5)==0)) {
        keys = (new tuple2<str *, str *>(2,key_text,key_text));
    }
    else {
        (__sys__::__ss_stderr)->write(const_6);
        return False;
    }
    len_key1 = len(keys->__getfirst__());
    len_key2 = len(keys->__getsecond__());
    handle_key = __OR(___bool((len_key1>0)), ___bool((len_key2>0)), 86);
    return True;
}

tuple2<str *, str *> *process_PE(str *read1, str *qual1, str *read2, str *qual2) {
    __iter<__ss_int> *__116;
    list<str *> *__109, *new_lseq;
    list<str *>::for_in_loop __112;
    tuple2<str *, __ss_int> *__159, *__162;
    __ss_bool __100, __101, __102, __103, __104, __105, __106, __107, __108, __113, __114, __124, __125, __126, __127, __137, __138, __139, __140, __141, __142, __143, __144, __145, __146, __153, __154, __88, __89, __90, __91, __92, __93, __94, __95, __96, __97, __98, __99, have_merged;
    __iter<str *> *__110;
    tuple2<str *, str *> *__123, *__136;
    list<__ss_int> *__115, *help_lqual2, *lqual1, *lqual2, *new_lqual, *rlqual2;
    list<double> *help_qual2, *pqual1, *pqual2, *rqual2;
    list<__ss_int>::for_in_loop __118;
    str *chimera, *help_rread2, *new_qual, *new_seq, *rread2;
    __ss_int __111, __117, __119, __120, __121, __122, __128, __129, __130, __131, __132, __133, __134, __135, __147, __148, __149, __150, __152, __156, __157, __158, __160, __161, clength, i, ipos, lread1, lread2, max_pos1, mlength, pos, start;
    double __151, __155, cval, cval1, cval2, hcval, max_value1;

    if (((len(read1)!=len(qual1)) or (len(read2)!=len(qual2)))) {
        (__sys__::__ss_stderr)->write(const_7);
        return (new tuple2<str *, str *>(3,const_2,const_2,const_2));
    }
    if ((handle_key and (len(read1)>0))) {
        if (((__eq(read1->__slice__(2, 0, len_key1, 0), keys->__getfirst__()) and __eq(read2->__slice__(2, 0, len_key2, 0), keys->__getsecond__())) or (options_allowMissing and (edits(read1->__slice__(2, 0, len_key1, 0), keys->__getfirst__())==1) and __eq(read2->__slice__(2, 0, len_key2, 0), keys->__getsecond__())) or (options_allowMissing and __eq(read1->__slice__(2, 0, len_key1, 0), keys->__getfirst__()) and (edits(read2->__slice__(2, 0, len_key2, 0), keys->__getsecond__())==1)))) {
            read1 = read1->__slice__(1, len_key1, 0, 0);
            qual1 = qual1->__slice__(1, len_key1, 0, 0);
            read2 = read2->__slice__(1, len_key2, 0, 0);
            qual2 = qual2->__slice__(1, len_key2, 0, 0);
        }
        else if ((options_allowMissing and __eq(read1->__slice__(2, 0, (len_key1-1), 0), (keys->__getfirst__())->__slice__(1, 1, 0, 0)) and __eq(read2->__slice__(2, 0, len_key2, 0), keys->__getsecond__()))) {
            read1 = read1->__slice__(1, (len_key1-1), 0, 0);
            qual1 = qual1->__slice__(1, (len_key1-1), 0, 0);
            read2 = read2->__slice__(1, len_key2, 0, 0);
            qual2 = qual2->__slice__(1, len_key2, 0, 0);
        }
        else if ((options_allowMissing and __eq(read1->__slice__(2, 0, len_key1, 0), keys->__getfirst__()) and __eq(read2->__slice__(2, 0, (len_key2-1), 0), (keys->__getsecond__())->__slice__(1, 1, 0, 0)))) {
            read1 = read1->__slice__(1, len_key1, 0, 0);
            qual1 = qual1->__slice__(1, len_key1, 0, 0);
            read2 = read2->__slice__(1, (len_key2-1), 0, 0);
            qual2 = qual2->__slice__(1, (len_key2-1), 0, 0);
        }
        else {
            return (new tuple2<str *, str *>(3,const_8,const_2,const_2));
        }
    }
    if ((len(adapter_chimeras)>0)) {
        lqual1 = convert_quality_logprob(qual1);
        pqual1 = convert_logprob_prob(lqual1);
        lread1 = len(read1);

        FOR_IN(chimera,adapter_chimeras,109,111,112)
            if ((quality_ident(read1, pqual1, chimera, default_0, maxadapter_comp)>cutoff_merge_trim)) {
                read1 = const_2;
                read2 = const_2;
                break;
            }
            else if ((quality_ident(read1, pqual1, chimera->__slice__(1, 1, 0, 0), default_0, maxadapter_comp)>cutoff_merge_trim)) {
                read1 = const_2;
                read2 = const_2;
                break;
            }
        END_FOR

        if ((__eq(read1, const_2) or __eq(read2, const_2))) {
            return (new tuple2<str *, str *>(3,const_9,const_2,const_2));
        }
    }
    lread1 = len(read1);
    lread2 = len(read2);
    lqual1 = convert_quality_logprob(qual1);
    pqual1 = convert_logprob_prob(lqual1);
    lqual2 = convert_quality_logprob(qual2);
    pqual2 = convert_logprob_prob(lqual2);
    rread2 = revcompl(read2);
    rqual2 = (new list<double>(pqual2));
    rqual2->reverse();
    rlqual2 = (new list<__ss_int>(lqual2));
    rlqual2->reverse();
    mlength = ___min(2, 0, lread1, lread2);
    clength = ___max(2, 0, lread1, lread2);
    have_merged = False;

    FOR_IN(start,(range((mlength+1)))->__slice__(4, 0, 0, (-1)),115,117,118)
        cval1 = 0.0;

        FAST_FOR(i,0,len(options_adapter_F),1,119,120)
            hcval = quality_ident(read1->__slice__(1, start, 0, 0), pqual1->__slice__(1, start, 0, 0), options_adapter_F->__getfast__(i), default_0, maxadapter_comp);
            if ((hcval>cval1)) {
                cval1 = hcval;
            }
        END_FOR

        cval2 = 0.0;

        FAST_FOR(i,0,len(options_adapter_S),1,121,122)
            hcval = quality_ident(read2->__slice__(1, start, 0, 0), pqual2->__slice__(1, start, 0, 0), options_adapter_S->__getfast__(i), default_0, maxadapter_comp);
            if ((hcval>cval2)) {
                cval2 = hcval;
            }
        END_FOR

        __123 = check_merge(read1->__slice__(2, 0, start, 0), pqual1->__slice__(2, 0, start, 0), lqual1->__slice__(2, 0, start, 0), rread2->__slice__(1, ___max(2, 0, 0, (lread2-start)), 0, 0), rqual2->__slice__(1, ___max(2, 0, 0, (lread2-start)), 0, 0), rlqual2->__slice__(1, ___max(2, 0, 0, (lread2-start)), 0, 0));
        new_seq = __123->__getfirst__();
        new_qual = __123->__getsecond__();
        if ((__ne(new_seq, const_2) and ((cval1>cutoff_merge_trim) or (cval2>cutoff_merge_trim)))) {
            have_merged = True;
            return (new tuple2<str *, str *>(3,const_2,new_seq,new_qual));
        }
    END_FOR

    if (__NOT(have_merged)) {

        FAST_FOR(start,(mlength+1),(clength+1),1,128,129)
            cval1 = 0.0;

            FAST_FOR(i,0,len(options_adapter_F),1,130,131)
                hcval = quality_ident(read1->__slice__(1, start, 0, 0), pqual1->__slice__(1, start, 0, 0), options_adapter_F->__getfast__(i), default_0, maxadapter_comp);
                if ((hcval>cval1)) {
                    cval1 = hcval;
                }
            END_FOR

            cval2 = 0.0;

            FAST_FOR(i,0,len(options_adapter_S),1,132,133)
                hcval = quality_ident(read2->__slice__(1, start, 0, 0), pqual2->__slice__(1, start, 0, 0), options_adapter_S->__getfast__(i), default_0, maxadapter_comp);
                if ((hcval>cval2)) {
                    cval2 = hcval;
                }
            END_FOR

            help_rread2 = const_2;
            help_lqual2 = (new list<__ss_int>());
            help_qual2 = (new list<double>());

            FAST_FOR(ipos,0,___max(2, 0, 0, (start-lread2)),1,134,135)
                help_rread2 = (help_rread2)->__iadd__(const_0);
                help_qual2->append(0.0);
                help_lqual2->append(0);
            END_FOR

            __136 = check_merge(read1->__slice__(2, 0, start, 0), pqual1->__slice__(2, 0, start, 0), lqual1->__slice__(2, 0, start, 0), (help_rread2)->__add__(rread2->__slice__(1, ___max(2, 0, 0, (lread2-start)), 0, 0)), (help_qual2)->__add__(rqual2->__slice__(1, ___max(2, 0, 0, (lread2-start)), 0, 0)), (help_lqual2)->__add__(rlqual2->__slice__(1, ___max(2, 0, 0, (lread2-start)), 0, 0)));
            new_seq = __136->__getfirst__();
            new_qual = __136->__getsecond__();
            if ((__ne(new_seq, const_2) and (((cval1>cutoff_merge_trim) or (cval2>cutoff_merge_trim)) or (__eq(read1->__slice__(1, start, 0, 0), read2->__slice__(1, start, 0, 0)) and (len(read1->__slice__(1, start, 0, 0))==0))))) {
                have_merged = True;
                if ((len(new_seq)<min_length)) {
                    return (new tuple2<str *, str *>(3,const_9,const_2,const_2));
                }
                else {
                    return (new tuple2<str *, str *>(3,const_2,new_seq,new_qual));
                }
                break;
            }
        END_FOR

    }
    if ((__NOT(have_merged) and options_mergeoverlap)) {
        __147 = (-1);
        __148 = (-1);
        max_value1 = __147;
        max_pos1 = __148;

        FAST_FOR(start,0,(lread1-options_min_overlap_seqs),1,149,150)
            if (((lread1-start)<=lread2)) {
                cval = quality_ident(read1->__slice__(1, start, 0, 0), pqual1->__slice__(1, start, 0, 0), rread2, rqual2, (-1));
                if ((cval>cutoff_merge_seqs_early)) {
                    __151 = cval;
                    __152 = start;
                    max_value1 = __151;
                    max_pos1 = __152;
                    break;
                }
                else if (((max_value1==((double)((-1)))) or (cval>max_value1))) {
                    __155 = cval;
                    __156 = start;
                    max_value1 = __155;
                    max_pos1 = __156;
                }
            }
        END_FOR

        if ((max_value1>cutoff_merge_seqs)) {
            if (options_onlyoverlap) {
                new_lseq = (new list<str *>(read1->__slice__(3, max_pos1, mlength, 0)));
                new_lqual = lqual1->__slice__(3, max_pos1, mlength, 0);

                FAST_FOR(pos,0,(mlength-max_pos1),1,157,158)
                    __159 = cons_base_prob(new_lseq->__getfast__(pos), read1->__getfast__((max_pos1+pos)), rqual2->__getfast__(pos), pqual1->__getfast__((max_pos1+pos)));
                    new_lseq->__setitem__(pos, __159->__getfirst__());
                    new_lqual->__setitem__(pos, __159->__getsecond__());
                END_FOR

            }
            else {
                new_lseq = (new list<str *>((read1->__slice__(2, 0, max_pos1, 0))->__add__(rread2)));
                new_lqual = (lqual1->__slice__(2, 0, max_pos1, 0))->__add__(rlqual2);

                FAST_FOR(pos,max_pos1,mlength,1,160,161)
                    __162 = cons_base_prob(new_lseq->__getfast__(pos), read1->__getfast__(pos), rqual2->__getfast__((pos-max_pos1)), pqual1->__getfast__(pos));
                    new_lseq->__setitem__(pos, __162->__getfirst__());
                    new_lqual->__setitem__(pos, __162->__getsecond__());
                END_FOR

            }
            return (new tuple2<str *, str *>(3,const_2,(const_2)->join(new_lseq),convert_logprob_quality(new_lqual)));
        }
    }
    return (new tuple2<str *, str *>(3,const_2,const_2,const_2));
}

tuple2<str *, str *> *overlap_reads(str *read1, str *qual1, str *read2, str *qual2) {
    list<str *> *new_lseq;
    double __169, __173, cval, max_value1;
    __ss_bool __163, __164, __171, __172;
    str *rread2;
    list<__ss_int> *lqual1, *lqual2, *new_lqual, *rlqual2;
    list<double> *pqual1, *pqual2, *rqual2;
    __ss_int __165, __166, __167, __168, __170, __174, __175, __176, lread1, lread2, max_pos1, mlength, pos, start;
    tuple2<str *, __ss_int> *__177;

    if (((len(read1)!=len(qual1)) or (len(read2)!=len(qual2)))) {
        (__sys__::__ss_stderr)->write(const_7);
        return (new tuple2<str *, str *>(3,const_2,const_2,const_2));
    }
    lread1 = len(read1);
    lread2 = len(read2);
    lqual1 = convert_quality_logprob(qual1);
    pqual1 = convert_logprob_prob(lqual1);
    lqual2 = convert_quality_logprob(qual2);
    pqual2 = convert_logprob_prob(lqual2);
    rread2 = revcompl(read2);
    rqual2 = (new list<double>(pqual2));
    rqual2->reverse();
    rlqual2 = (new list<__ss_int>(lqual2));
    rlqual2->reverse();
    mlength = ___min(2, 0, lread1, lread2);
    __165 = (-1);
    __166 = (-1);
    max_value1 = __165;
    max_pos1 = __166;

    FAST_FOR(start,0,(mlength-options_min_overlap_seqs),1,167,168)
        if (((lread2-start)<=lread1)) {
            cval = quality_ident(read1, pqual1, rread2->__slice__(1, start, 0, 0), rqual2->__slice__(1, start, 0, 0), (-1));
            if ((cval>cutoff_merge_seqs_early)) {
                __169 = cval;
                __170 = start;
                max_value1 = __169;
                max_pos1 = __170;
                break;
            }
            else if (((max_value1==((double)((-1)))) or (cval>max_value1))) {
                __173 = cval;
                __174 = start;
                max_value1 = __173;
                max_pos1 = __174;
            }
        }
    END_FOR

    if ((max_value1>cutoff_merge_seqs)) {
        new_lseq = (new list<str *>(rread2->__slice__(3, max_pos1, mlength, 0)));
        new_lqual = rlqual2->__slice__(3, max_pos1, mlength, 0);

        FAST_FOR(pos,0,(mlength-max_pos1),1,175,176)
            __177 = cons_base_prob(new_lseq->__getfast__(pos), read1->__getfast__(pos), rqual2->__getfast__((max_pos1+pos)), pqual1->__getfast__(pos));
            new_lseq->__setitem__(pos, __177->__getfirst__());
            new_lqual->__setitem__(pos, __177->__getsecond__());
        END_FOR

        return (new tuple2<str *, str *>(3,const_2,(const_2)->join(new_lseq),convert_logprob_quality(new_lqual)));
    }
    else {
        return (new tuple2<str *, str *>(3,const_2,const_2,const_2));
    }
    return 0;
}

tuple2<str *, str *> *consensus_reads(str *read1, str *qual1, str *read2, str *qual2) {
    list<str *> *new_lseq;
    tuple2<str *, __ss_int> *__183;
    __ss_bool __178, __179, __180;
    list<__ss_int> *lqual1, *lqual2, *new_lqual;
    list<double> *pqual1, *pqual2;
    __ss_int __181, __182, pos;
    double cval;

    if (((len(read1)!=len(qual1)) or (len(read2)!=len(qual2)) or (len(read1)!=len(read2)))) {
        (__sys__::__ss_stderr)->write(const_7);
        return (new tuple2<str *, str *>(3,const_2,const_2,const_2));
    }
    lqual1 = convert_quality_logprob(qual1);
    pqual1 = convert_logprob_prob(lqual1);
    lqual2 = convert_quality_logprob(qual2);
    pqual2 = convert_logprob_prob(lqual2);
    cval = quality_ident(read1, pqual1, read2, pqual2, (-1));
    if ((cval>0.75)) {
        new_lseq = (new list<str *>(read2));
        new_lqual = lqual2;

        FAST_FOR(pos,0,len(read1),1,181,182)
            __183 = cons_base_prob(new_lseq->__getfast__(pos), read1->__getfast__(pos), pqual2->__getfast__(pos), pqual1->__getfast__(pos));
            new_lseq->__setitem__(pos, __183->__getfirst__());
            new_lqual->__setitem__(pos, __183->__getsecond__());
        END_FOR

        return (new tuple2<str *, str *>(3,const_2,(const_2)->join(new_lseq),convert_logprob_quality(new_lqual)));
    }
    else {
        return (new tuple2<str *, str *>(3,const_2,const_2,const_2));
    }
    return 0;
}

tuple2<str *, str *> *process_SR(str *read1, str *qual1) {
    list<str *> *__190;
    list<str *>::for_in_loop __193;
    double __200, __209, cval, hcval, max_value;
    __ss_bool __184, __185, __186, __187, __188, __189, __202, __203, __204, __205, __206, __207, __208, __211, __212, key_handled;
    str *chimera;
    __iter<str *> *__191;
    list<double> *pqual1;
    list<__ss_int> *lqual1;
    __ss_int __192, __194, __195, __196, __197, __198, __199, __201, __210, adapter_pos, i, lread1, max_pos, start;

    if ((len(read1)!=len(qual1))) {
        (__sys__::__ss_stderr)->write(const_7);
        return (new tuple2<str *, str *>(3,const_2,const_2,const_2));
    }
    if (handle_key) {
        key_handled = False;
        if ((__eq(read1->__slice__(2, 0, len_key1, 0), keys->__getfirst__()) or (options_allowMissing and (edits(read1->__slice__(2, 0, len_key1, 0), keys->__getfirst__())==1)))) {
            read1 = read1->__slice__(1, len_key1, 0, 0);
            qual1 = qual1->__slice__(1, len_key1, 0, 0);
            key_handled = True;
        }
        else if ((options_allowMissing and __eq(read1->__slice__(2, 0, (len_key1-1), 0), (keys->__getfirst__())->__slice__(1, 1, 0, 0)))) {
            read1 = read1->__slice__(1, (len_key1-1), 0, 0);
            qual1 = qual1->__slice__(1, (len_key1-1), 0, 0);
            key_handled = True;
        }
        if (__NOT(key_handled)) {
            return (new tuple2<str *, str *>(3,const_8,const_2,const_2));
        }
    }
    lqual1 = convert_quality_logprob(qual1);
    pqual1 = convert_logprob_prob(lqual1);
    lread1 = len(read1);
    if ((len(adapter_chimeras)>0)) {

        FOR_IN(chimera,adapter_chimeras,190,192,193)
            if ((quality_ident(read1, pqual1, chimera, default_0, maxadapter_comp)>cutoff_merge_trim)) {
                read1 = const_2;
                break;
            }
            else if ((quality_ident(read1, pqual1, chimera->__slice__(1, 1, 0, 0), default_0, (maxadapter_comp-1))>cutoff_merge_trim)) {
                read1 = const_2;
                break;
            }
        END_FOR

        if (__eq(read1, const_2)) {
            return (new tuple2<str *, str *>(3,const_9,const_2,const_2));
        }
    }
    adapter_pos = (-2);
    __194 = (-1);
    __195 = (-1);
    max_value = __194;
    max_pos = __195;

    FAST_FOR(start,0,lread1,1,196,197)
        cval = 0.0;

        FAST_FOR(i,0,len(options_adapter_F),1,198,199)
            hcval = quality_ident(read1->__slice__(1, start, 0, 0), pqual1->__slice__(1, start, 0, 0), options_adapter_F->__getfast__(i), default_0, maxadapter_comp);
            if ((hcval>cval)) {
                cval = hcval;
            }
        END_FOR

        if ((cval>cutoff_merge_seqs_early)) {
            __200 = cval;
            __201 = start;
            max_value = __200;
            max_pos = __201;
            break;
        }
        else if ((((cval>cutoff_merge_trim) and (cval>max_value) and ((lread1-start)>min_length)) or (((lread1-start)<=min_length) and (max_value==((double)((-1))))))) {
            __209 = cval;
            __210 = start;
            max_value = __209;
            max_pos = __210;
        }
    END_FOR

    if (((max_value>cutoff_merge_trim) and ((lread1-max_pos)>=options_trimCutoff))) {
        adapter_pos = max_pos;
        read1 = read1->__slice__(2, 0, adapter_pos, 0);
        if ((len(read1)<min_length)) {
            return (new tuple2<str *, str *>(3,const_9,const_2,const_2));
        }
        else {
            return (new tuple2<str *, str *>(3,const_2,read1,qual1->__slice__(2, 0, adapter_pos, 0)));
        }
    }
    return (new tuple2<str *, str *>(3,const_2,const_2,const_2));
}

void __init() {
    const_0 = __char_cache[73];;
    const_1 = new str("Quality_ident: Call with invalid arguments!\n");
    const_2 = new str("");
    const_3 = new str("ACGT");
    const_4 = __char_cache[78];;
    const_5 = __char_cache[44];;
    const_6 = new str("Unexpected number of keys specified.\n");
    const_7 = new str("Error: Length of quality score and read strings do not match!\n");
    const_8 = __char_cache[75];;
    const_9 = __char_cache[68];;
    const_10 = new str("TGCA");
    const_11 = new str("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG");
    const_12 = new str("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT");
    const_13 = new str("ACACTCTTTCCCTACACGTCTGAACTCCAG,ACACTCTTTCCCACACGTCTGAACTCCAGT,ACACTCTTTCCCTACACACGTCTGAACTCC,CTCTTTCCCTACACGTCTGAACTCCAGTCA,GAAGAGCACACGTCTGAACTCCAGTCACII,GAGCACACGTCTGAACTCCAGTCACIIIII,GATCGGAAGAGCACACGTCTGAACTCCAGT,AGATCGGAAGAGCACACGTCTGAACTCCAG,AGAGCACACGTCTGAACTCCAGTCACIIII,ACACGTCTGAACTCCAGTCACIIIIIIIAT,GTGCACACGTCTGAACTCCAGTCACIIIII,AGCACACGTCTGAACTCCAGTCACIIIIII,CGTATGCCGTCTTCTGCTTGAAAAAAAAAA");
    const_14 = new str("__main__");
    const_15 = new str("Read cosensus");
    const_16 = new str("ACGTACGTACGT");
    const_17 = new str("000000000000");
    const_18 = new str("ACGTACGTCCGT");
    const_19 = new str("SR");
    const_20 = new str("PE fail");
    const_21 = new str("TGCATGCATGCA");
    const_22 = new str("PE no adapter");
    const_23 = new str("PE adapter");
    const_24 = new str("ATAAACATATGGCAAACATGGTTCTAGATCGGAAGAGCACACGT");
    const_25 = new str("00000000000000000000000000000000000000000000");
    const_26 = new str("AGAACCATGTTTGCCATATGTTTATAGATCGGAAGAGCGTCGTG");
    const_27 = new str("PE overlap consensus");
    const_28 = new str("ATAAACATATGGCAAACATGGTTCTAAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAATCACTAACACATTCATTTTCCAGACACCCTTCACACTACCGTCGGATCGTGCGTGTACCTCTGAATCTCGTATGCCGTCTTCT");
    const_29 = new str("55???BBBDDDDDDDDFFFFF>EEHBF?EFFFFG@FHHHHHHHHHHHH@GHHHHHHHHHHHHHHF?FFHEGDG=EHG?FFGCFGBGHHHFGCFCFHHHHHBD.?C--@=<+<CECE=;46=@@+=,=,5=5,AB,,55A;*5,C=:;18?##");
    const_30 = new str("TCTGGAAAATGAATGTGTTAGTGATTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTAGAACCATGTTTGCCATATGTTTATCTTCAGCTTCCCGATTACGGATCTCGTATGTGTAGATCTCGGTGGTCGCCG");
    const_31 = new str("+?BDDDDDFFFCCBFEFHFFF>EF@GGHHHHFHFHHHFFHHHHHHFFHHHHHGHHHHHHHHHHHHHHH-A-CFF,C//AECFD?EEFDFGH/CAFHF/ACFHDFFEEEHDH<DDD)7?;?+@7@-7,66B?D?,+@@@######");
    const_32 = new str("Overlap read consensus");
    const_33 = new str("ATTGCGTGAACCGAAGCTCATCAAGATCTG");
    const_34 = new str("GTTGATCCGGTCCTAGGCAGTGAAGATCTC");
    const_35 = new str("TTTACGGCTCATTGCGTGAACCGAAGCTCATCAAGATCTGGCCTCGGCGGCCAAGCTTAGTCGCCTATACGGTGATGGGTGCATNNTTCATNNNNATNNNCTCCCCGTTTAATCCCATATCTCGTATGCCGTCTTCTGCTTGAAAAAAAA");
    const_36 = new str("=======++5<5@9@@>CCEEE>CCCCCEEFFGFFFFFFDFFEFEEEECCDDCE+ACEDD5CEDEEECDE@D=9E+CDD@;@DE##113D9####00###21*28@@2<E?E=E=EEE</;(6?EEE<E<<E;6<EE#############");
    const_37 = new str("AGCCGTAAAAGTTGATCCGGTCCTAGGCAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAAAAAAAACAAGGAAAATAAAAAAAAAAACCCAGCCAACTAACAAAAAACAAACAAAACAAAAACGAAACCAAAACCAAAAAATAC");
    const_38 = new str("????????BDB<9<BBFFBCCC7ECCACEECFGFCDBCBFFFEEDC?>7C7>EEFGH=CDEGGHHHDE##################################################################################");

    __name__ = new str("MergeTrimReads");

    table = __string__::maketrans(const_10, const_3);
    offset = 33;
    cutoff_merge_trim = 0.8;
    cutoff_merge_seqs_early = 0.95;
    cutoff_merge_seqs = 0.9;
    options_min_overlap_seqs = 10;
    max_prob_same_channel = 0.5;
    max_prob_N = 0.25;
    maxadapter_comp = 30;
    min_length = 5;
    options_adapter_F = (new list<str *>(1,const_11));
    options_adapter_S = (new list<str *>(1,const_12));
    options_adapter_chimera = const_13;
    adapter_chimeras = options_adapter_chimera->split(const_5);
    options_trimCutoff = 1;
    keys = (new tuple2<str *, str *>(2,const_2,const_2));
    len_key1 = len(keys->__getfirst__());
    len_key2 = len(keys->__getsecond__());
    handle_key = __OR(___bool((len_key1>0)), ___bool((len_key2>0)), 3);
    options_allowMissing = False;
    options_mergeoverlap = False;
    options_onlyoverlap = False;
    default_0 = (new list<double>());
    default_1 = False;
    default_2 = False;
    default_3 = False;
    default_4 = const_2;
    default_5 = const_2;
    default_6 = const_2;
    if (__eq(__name__, const_14)) {
        set_adapter_sequences((const_5)->join(options_adapter_F), (const_5)->join(options_adapter_S), options_adapter_chimera, maxadapter_comp);
        set_options(options_trimCutoff, options_allowMissing, True, False, 10);
        set_keys(const_5);
        print2(NULL,0,2, const_15, consensus_reads(const_16, const_17, const_18, const_17));
        print2(NULL,0,2, const_19, process_SR(const_16, const_17));
        print2(NULL,0,2, const_20, process_PE(const_16, const_17, const_21, const_17));
        print2(NULL,0,2, const_22, process_PE(const_16, const_17, const_16, const_17));
        print2(NULL,0,2, const_23, process_PE(const_24, const_25, const_26, const_25));
        print2(NULL,0,2, const_27, process_PE(const_28, const_29, const_30, const_31));
        print2(NULL,0,2, const_32, overlap_reads(const_28, const_29, const_30, const_31));
        print2(NULL,0,0);
        set_adapter_sequences(const_33, const_34, const_2, maxadapter_comp);
        print2(NULL,0,0);
        print2(NULL,0,1, process_SR(const_35, const_36));
        set_adapter_sequences(const_34, const_33, const_2, maxadapter_comp);
        print2(NULL,0,1, process_SR(const_37, const_38));
        print2(NULL,0,0);
        print2(NULL,0,1, process_PE(const_35, const_36, const_37, const_38));
    }
}

} // module namespace

/* extension module glue */

extern "C" {
#include <Python.h>
#include "os/path.hpp"
#include "os/__init__.hpp"
#include "string.hpp"
#include "sys.hpp"
#include "math.hpp"
#include "stat.hpp"
#include "MergeTrimReads.hpp"
#include <structmember.h>
#include "os/path.hpp"
#include "os/__init__.hpp"
#include "string.hpp"
#include "sys.hpp"
#include "math.hpp"
#include "stat.hpp"
#include "MergeTrimReads.hpp"

PyObject *__ss_mod_MergeTrimReads;

namespace __MergeTrimReads__ { /* XXX */
PyObject *Global_MergeTrimReads_quality_ident(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        str *arg_0 = __ss_arg<str *>("seq1", 0, 0, 0, args, kwargs);
        list<double> *arg_1 = __ss_arg<list<double> *>("qual1", 1, 0, 0, args, kwargs);
        str *arg_2 = __ss_arg<str *>("seq2", 2, 0, 0, args, kwargs);
        list<double> *arg_3 = __ss_arg<list<double> *>("qual2", 3, 1, __MergeTrimReads__::default_0, args, kwargs);
        __ss_int arg_4 = __ss_arg<__ss_int >("maxcomp", 4, 1, (-1), args, kwargs);

        return __to_py(__MergeTrimReads__::quality_ident(arg_0, arg_1, arg_2, arg_3, arg_4));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->message)?(e->message->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_MergeTrimReads_set_adapter_sequences(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        str *arg_0 = __ss_arg<str *>("forward", 0, 1, __MergeTrimReads__::default_4, args, kwargs);
        str *arg_1 = __ss_arg<str *>("reverse", 1, 1, __MergeTrimReads__::default_5, args, kwargs);
        str *arg_2 = __ss_arg<str *>("chimera", 2, 1, __MergeTrimReads__::default_6, args, kwargs);
        __ss_int arg_3 = __ss_arg<__ss_int >("max_comp", 3, 1, 30, args, kwargs);

        return __to_py(__MergeTrimReads__::set_adapter_sequences(arg_0, arg_1, arg_2, arg_3));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->message)?(e->message->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_MergeTrimReads_set_options(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        __ss_int arg_0 = __ss_arg<__ss_int >("trimcutoff", 0, 1, 1, args, kwargs);
        __ss_bool arg_1 = __ss_arg<__ss_bool >("allowMissing", 1, 1, __MergeTrimReads__::default_1, args, kwargs);
        __ss_bool arg_2 = __ss_arg<__ss_bool >("mergeoverlap", 2, 1, __MergeTrimReads__::default_2, args, kwargs);
        __ss_bool arg_3 = __ss_arg<__ss_bool >("onlyoverlap", 3, 1, __MergeTrimReads__::default_3, args, kwargs);
        __ss_int arg_4 = __ss_arg<__ss_int >("min_overlap_seqs", 4, 1, 10, args, kwargs);

        return __to_py(__MergeTrimReads__::set_options(arg_0, arg_1, arg_2, arg_3, arg_4));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->message)?(e->message->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_MergeTrimReads_revcompl(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        str *arg_0 = __ss_arg<str *>("seq", 0, 0, 0, args, kwargs);

        return __to_py(__MergeTrimReads__::revcompl(arg_0));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->message)?(e->message->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_MergeTrimReads_process_SR(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        str *arg_0 = __ss_arg<str *>("read1", 0, 0, 0, args, kwargs);
        str *arg_1 = __ss_arg<str *>("qual1", 1, 0, 0, args, kwargs);

        return __to_py(__MergeTrimReads__::process_SR(arg_0, arg_1));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->message)?(e->message->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_MergeTrimReads_convert_quality_logprob(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        str *arg_0 = __ss_arg<str *>("qualstring", 0, 0, 0, args, kwargs);

        return __to_py(__MergeTrimReads__::convert_quality_logprob(arg_0));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->message)?(e->message->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_MergeTrimReads_edits(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        str *arg_0 = __ss_arg<str *>("seq1", 0, 0, 0, args, kwargs);
        str *arg_1 = __ss_arg<str *>("seq2", 1, 0, 0, args, kwargs);

        return __to_py(__MergeTrimReads__::edits(arg_0, arg_1));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->message)?(e->message->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_MergeTrimReads_cons_base_prob(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        str *arg_0 = __ss_arg<str *>("base1", 0, 0, 0, args, kwargs);
        str *arg_1 = __ss_arg<str *>("base2", 1, 0, 0, args, kwargs);
        double arg_2 = __ss_arg<double >("prob1", 2, 0, 0, args, kwargs);
        double arg_3 = __ss_arg<double >("prob2", 3, 0, 0, args, kwargs);

        return __to_py(__MergeTrimReads__::cons_base_prob(arg_0, arg_1, arg_2, arg_3));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->message)?(e->message->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_MergeTrimReads_overlap_reads(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        str *arg_0 = __ss_arg<str *>("read1", 0, 0, 0, args, kwargs);
        str *arg_1 = __ss_arg<str *>("qual1", 1, 0, 0, args, kwargs);
        str *arg_2 = __ss_arg<str *>("read2", 2, 0, 0, args, kwargs);
        str *arg_3 = __ss_arg<str *>("qual2", 3, 0, 0, args, kwargs);

        return __to_py(__MergeTrimReads__::overlap_reads(arg_0, arg_1, arg_2, arg_3));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->message)?(e->message->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_MergeTrimReads_check_merge(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        str *arg_0 = __ss_arg<str *>("read1", 0, 0, 0, args, kwargs);
        list<double> *arg_1 = __ss_arg<list<double> *>("qual1", 1, 0, 0, args, kwargs);
        list<__ss_int> *arg_2 = __ss_arg<list<__ss_int> *>("pqual1", 2, 0, 0, args, kwargs);
        str *arg_3 = __ss_arg<str *>("read2", 3, 0, 0, args, kwargs);
        list<double> *arg_4 = __ss_arg<list<double> *>("qual2", 4, 0, 0, args, kwargs);
        list<__ss_int> *arg_5 = __ss_arg<list<__ss_int> *>("pqual2", 5, 0, 0, args, kwargs);

        return __to_py(__MergeTrimReads__::check_merge(arg_0, arg_1, arg_2, arg_3, arg_4, arg_5));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->message)?(e->message->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_MergeTrimReads_convert_logprob_prob(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        list<__ss_int> *arg_0 = __ss_arg<list<__ss_int> *>("lprobs", 0, 0, 0, args, kwargs);

        return __to_py(__MergeTrimReads__::convert_logprob_prob(arg_0));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->message)?(e->message->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_MergeTrimReads_convert_logprob_quality(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        list<__ss_int> *arg_0 = __ss_arg<list<__ss_int> *>("probs", 0, 0, 0, args, kwargs);

        return __to_py(__MergeTrimReads__::convert_logprob_quality(arg_0));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->message)?(e->message->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_MergeTrimReads_consensus_reads(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        str *arg_0 = __ss_arg<str *>("read1", 0, 0, 0, args, kwargs);
        str *arg_1 = __ss_arg<str *>("qual1", 1, 0, 0, args, kwargs);
        str *arg_2 = __ss_arg<str *>("read2", 2, 0, 0, args, kwargs);
        str *arg_3 = __ss_arg<str *>("qual2", 3, 0, 0, args, kwargs);

        return __to_py(__MergeTrimReads__::consensus_reads(arg_0, arg_1, arg_2, arg_3));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->message)?(e->message->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_MergeTrimReads_process_PE(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        str *arg_0 = __ss_arg<str *>("read1", 0, 0, 0, args, kwargs);
        str *arg_1 = __ss_arg<str *>("qual1", 1, 0, 0, args, kwargs);
        str *arg_2 = __ss_arg<str *>("read2", 2, 0, 0, args, kwargs);
        str *arg_3 = __ss_arg<str *>("qual2", 3, 0, 0, args, kwargs);

        return __to_py(__MergeTrimReads__::process_PE(arg_0, arg_1, arg_2, arg_3));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->message)?(e->message->unit.c_str()):""));
        return 0;
    }
}

PyObject *Global_MergeTrimReads_set_keys(PyObject *self, PyObject *args, PyObject *kwargs) {
    try {
        str *arg_0 = __ss_arg<str *>("key_text", 0, 0, 0, args, kwargs);

        return __to_py(__MergeTrimReads__::set_keys(arg_0));

    } catch (Exception *e) {
        PyErr_SetString(__to_py(e), ((e->message)?(e->message->unit.c_str()):""));
        return 0;
    }
}

static PyNumberMethods Global_MergeTrimReads_as_number = {
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
};

static PyMethodDef Global_MergeTrimReadsMethods[] = {
    {(char *)"__newobj__", (PyCFunction)__ss__newobj__, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"quality_ident", (PyCFunction)Global_MergeTrimReads_quality_ident, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"set_adapter_sequences", (PyCFunction)Global_MergeTrimReads_set_adapter_sequences, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"set_options", (PyCFunction)Global_MergeTrimReads_set_options, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"revcompl", (PyCFunction)Global_MergeTrimReads_revcompl, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"process_SR", (PyCFunction)Global_MergeTrimReads_process_SR, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"convert_quality_logprob", (PyCFunction)Global_MergeTrimReads_convert_quality_logprob, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"edits", (PyCFunction)Global_MergeTrimReads_edits, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"cons_base_prob", (PyCFunction)Global_MergeTrimReads_cons_base_prob, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"overlap_reads", (PyCFunction)Global_MergeTrimReads_overlap_reads, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"check_merge", (PyCFunction)Global_MergeTrimReads_check_merge, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"convert_logprob_prob", (PyCFunction)Global_MergeTrimReads_convert_logprob_prob, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"convert_logprob_quality", (PyCFunction)Global_MergeTrimReads_convert_logprob_quality, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"consensus_reads", (PyCFunction)Global_MergeTrimReads_consensus_reads, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"process_PE", (PyCFunction)Global_MergeTrimReads_process_PE, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {(char *)"set_keys", (PyCFunction)Global_MergeTrimReads_set_keys, METH_VARARGS | METH_KEYWORDS, (char *)""},
    {NULL}
};

PyMODINIT_FUNC initMergeTrimReads(void) {
    __shedskin__::__init();
    __sys__::__init(0, 0);
    __stat__::__init();
    __os__::__path__::__init();
    __os__::__init();
    __math__::__init();
    __string__::__init();
    __MergeTrimReads__::__init();

    __ss_mod_MergeTrimReads = Py_InitModule((char *)"MergeTrimReads", Global_MergeTrimReadsMethods);
    if(!__ss_mod_MergeTrimReads)
        return;


    addMergeTrimReads();

}

PyMODINIT_FUNC addMergeTrimReads(void) {
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"cutoff_merge_trim", __to_py(__MergeTrimReads__::cutoff_merge_trim));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"len_key2", __to_py(__MergeTrimReads__::len_key2));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"len_key1", __to_py(__MergeTrimReads__::len_key1));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"options_allowMissing", __to_py(__MergeTrimReads__::options_allowMissing));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"table", __to_py(__MergeTrimReads__::table));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"adapter_chimeras", __to_py(__MergeTrimReads__::adapter_chimeras));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"cutoff_merge_seqs_early", __to_py(__MergeTrimReads__::cutoff_merge_seqs_early));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"max_prob_N", __to_py(__MergeTrimReads__::max_prob_N));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"options_trimCutoff", __to_py(__MergeTrimReads__::options_trimCutoff));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"maxadapter_comp", __to_py(__MergeTrimReads__::maxadapter_comp));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"options_adapter_chimera", __to_py(__MergeTrimReads__::options_adapter_chimera));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"min_length", __to_py(__MergeTrimReads__::min_length));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"keys", __to_py(__MergeTrimReads__::keys));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"options_min_overlap_seqs", __to_py(__MergeTrimReads__::options_min_overlap_seqs));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"offset", __to_py(__MergeTrimReads__::offset));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"cutoff_merge_seqs", __to_py(__MergeTrimReads__::cutoff_merge_seqs));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"max_prob_same_channel", __to_py(__MergeTrimReads__::max_prob_same_channel));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"options_adapter_S", __to_py(__MergeTrimReads__::options_adapter_S));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"options_mergeoverlap", __to_py(__MergeTrimReads__::options_mergeoverlap));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"options_onlyoverlap", __to_py(__MergeTrimReads__::options_onlyoverlap));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"options_adapter_F", __to_py(__MergeTrimReads__::options_adapter_F));
    PyModule_AddObject(__ss_mod_MergeTrimReads, (char *)"handle_key", __to_py(__MergeTrimReads__::handle_key));

}

} // namespace __MergeTrimReads__

} // extern "C"
int main(int __ss_argc, char **__ss_argv) {
    __shedskin__::__init();
    __sys__::__init(0, 0);
    __stat__::__init();
    __os__::__path__::__init();
    __os__::__init();
    __math__::__init();
    __string__::__init();
    __shedskin__::__start(__MergeTrimReads__::__init);
}
