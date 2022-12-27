/* MIT License
 *
 * Copyright (c) 2018 Sam Kovaka <skovaka@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef _INCL_BWAFMI
#define _INCL_BWAFMI

#include <string>
#include <climits>
#include <utility>
#include <cstring>
#include <bwa/bwa.h>
#include <bwa/utils.h>
#include <zlib.h>
#include <algorithm>
#include "util.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#endif

//From submods/bwa/bwtindex.c
#define BWA_BLOCK_SIZE 10000000

using Range = std::pair<i64, i64>;
using GenomeChrs = std::vector< std::vector<std::string> >;

struct PanKmerBitvec {
    std::vector<i64> vec;
    //GenomeChrs chrs;
    size_t W,H;

    PanKmerBitvec(size_t seq_len, size_t genome_count, size_t k) : 
            W(seq_len - k + 1), H(genome_count) {
        vec.resize(W*H);
    }

    void set_genome(size_t i, size_t j) {
        vec[(i * W) + (j % W)] += 1;
    }
};

class RefIndex {
    public:

    RefIndex() :
        index_(NULL),
        bns_(NULL),
        pacseq_(NULL),
        loaded_(false),
        size_(0) {}

    RefIndex(const std::string &prefix, bool pacseq=false, bool bwt=true) : 
            index_(NULL),
            bns_(NULL),
            pacseq_(NULL),
            loaded_(false),
            size_(0) {
        if (!prefix.empty()) {
            bns_ = bns_restore(prefix.c_str());
            size_ = 2 * (bns_->l_pac);
            if (bwt) load_index(prefix);
            if (pacseq) load_pacseq();
        }
    }

    void load_index(const std::string &prefix) {
        std::string bwt_fname = prefix + ".bwt",
                    sa_fname = prefix + ".sa";

        if (bns_ == NULL) {
            bns_ = bns_restore(prefix.c_str());
        }

        index_ = bwt_restore_bwt(bwt_fname.c_str());
        bwt_restore_sa(sa_fname.c_str(), index_);

        //for (KmerType k = 0; k < kmer_ranges_.size(); k++) {
        //    Range r = get_base_range(ModelType::kmer_head(k));
        //    for (u8 i = 1; i < K; i++) {
        //        r = get_neighbor(r, ModelType::kmer_base(k, i));
        //    }
        //    kmer_ranges_[k] = r;
        //}

        loaded_ = true;
    }

    static void create(const std::string &fasta_fname, 
                       const std::string &prefix = "",
                       bool no_bwt=false) {

        std::string prefix_auto = prefix.empty() ? fasta_fname : prefix;

        if (no_bwt) {
            gzFile fp = xzopen(fasta_fname.c_str(), "r");
            bns_fasta2bntseq(fp, prefix_auto.c_str(), 0);
            err_gzclose(fp);

        } else {
            bwa_idx_build(fasta_fname.c_str(), 
                          prefix.c_str(), 
                          BWTALGO_AUTO,
                          BWA_BLOCK_SIZE);
        }
    }

    bool bwt_loaded() {
        return loaded_;
    }

    void load_pacseq() {
        if (!pacseq_loaded()) {
            //Copied from bwa/bwase.c
            pacseq_ = (u8*) calloc(bns_->l_pac/4+1, 1);
            err_fread_noeof(pacseq_, 1, bns_->l_pac/4+1, bns_->fp_pac);
        }   
    }

    void destroy() {
        if (index_ != NULL) { 
            bwt_destroy(index_);
        }
        if (bns_ != NULL) { 
            bns_destroy(bns_);
        }
        if (pacseq_loaded()) {
            free(pacseq_);
        }
    }

    Range get_neighbor(Range r1, u8 base) const {
        u64 os, oe;
        bwt_2occ(index_, r1.first - 1, r1.second, base, &os, &oe);
        return Range(index_->L2[base] + os + 1, index_->L2[base] + oe);
    }

    Range get_base_range(u8 base) const {
        return Range(index_->L2[base], index_->L2[base+1]);
    }

    i64 fm_to_mref(i64 fm) {
        return size() - fm_to_pac(fm);
    }

    i64 fm_to_pac(i64 fm) const {
        return bwt_sa(index_, fm);
    }

    size_t size() const {
        //return index_->seq_len;
        return size_;
    }

    i32 get_ref_id(const std::string &ref_name) {
        for (int i = 0; i < bns_->n_seqs; i++) {
            if (strcmp(bns_->anns[i].name, ref_name.c_str()) == 0) {
                return i;
            }
        }
        return -1;
    }

    i32 mref_to_ref_id(i64 mref) {
        return pac_to_ref_id(mref_to_pac(mref));
    }

    i32 pac_to_ref_id(i64 pac) {
        return bns_pos2rid(bns_, pac);
    }

    bool is_mref_flipped(i64 i) const {
        return i >= static_cast<i64>(size() / 2);
    }

    bool is_mref_fwd(i64 i, bool is_rna) const {
        return is_mref_flipped(i) == is_rna;
    }

    i64 mref_to_pac(i64 mref) {
        if (is_mref_flipped(mref)) {
            return size() - mref - 1;
        }
        return mref;
    }

    std::string get_ref_name(u32 rid) {
        return bns_->anns[rid].name;
    }


    i64 get_ref_len(u32 rid) const {
        return bns_->anns[rid].len;
    }

    i64 get_sa_loc(const std::string &name, i64 coord) {
        for (int i = 0; i < bns_->n_seqs; i++) {
            if (strcmp(bns_->anns[i].name, name.c_str()) == 0) {
                return bns_->anns[i].offset + coord;
            }
        }
        return 0;
    }


    std::vector< std::pair<std::string, i64> > get_seqs() const {
        std::vector< std::pair<std::string, i64> > seqs;

        for (i32 i = 0; i < bns_->n_seqs; i++) {
            bntann1_t ann = bns_->anns[i];
            std::string name = std::string(ann.name);
            seqs.push_back( std::pair<std::string, i64>(name, ann.len) );
        }

        return seqs;
    }

    i64 ref_to_pac(std::string name, i64 coord) {
        i32 i;
        for (i = 0; i < bns_->n_seqs; i++) {
            if (strcmp(name.c_str(), bns_->anns[i].name) == 0)
                return bns_->anns[i].offset + coord;
        }
        return INT_MAX;
    } 

    //TODO overload same function below
    i64 get_pac_shift(const std::string &ref_name) {
        return get_ref_shift(get_ref_id(ref_name));
    }

    i64 get_ref_shift(i32 ref_id) {
        return bns_->anns[ref_id].offset;
    }

    bool pacseq_loaded() const {
        return pacseq_ != NULL;
    }

    std::pair<i64, i64> ref_to_mref(const std::string &ref_name, i64 st, i64 en, bool is_fwd, bool is_rna) {
        auto shift = get_pac_shift(ref_name);

        i64 pac_st = shift+st, 
            pac_en = shift+en;

        auto flip = is_fwd == is_rna;

        if (!flip) return {pac_st, pac_en};
        return {size() - pac_en, size() - pac_st};
    }

    i64 ref_to_mref(i32 rid, i64 ref_coord, bool is_fwd, bool is_rna) {
        auto shift = bns_->anns[rid].offset;
        auto pac_coord = shift+ref_coord;
        auto flip = is_fwd == is_rna;

        if (!flip) return pac_coord;
        return size() - pac_coord - 1;
    }

    i64 mref_to_ref(i64 mref) {
        i64 pac;
        if (is_mref_flipped(mref)) {
            pac = size() - mref - 1;
        } else {
            pac = mref;
        }
        i32 rid = bns_pos2rid(bns_, pac);
        return pac - bns_->anns[rid].offset;
    }

    i64 pac_to_ref(i64 pac) {
        i32 rid = bns_pos2rid(bns_, pac);
        return pac - bns_->anns[rid].offset;
    }

    u8 get_base(i64 pac, bool comp=false) {
        if (pac < 0 || pac > static_cast<i64>(size() / 2)) { //TODO better size definition
            throw std::out_of_range("Base out of range: " + std::to_string(pac));
        }
        size_t i = pac >> 2,
               shift = ((pac & 3) ^ 3) << 1;
        return ((pacseq_[i] >> shift) & 3) ^ (comp | (comp<<1));
    }

    PanKmerBitvec pan_kmer_bitvec(std::vector<u8> seq, size_t k, GenomeChrs chr_groups) {
        //std::vector< std::vector<i32> > rids;
        std::unordered_map<int, int> rid_groups;
        for (size_t i = 0; i < chr_groups.size(); i++) {
            //rids.push_back({});
            for (auto &chr : chr_groups[i]) {
                rid_groups[get_ref_id(chr)] = i;
            }
        }

        PanKmerBitvec ret(seq.size(), chr_groups.size(), k);

        for (size_t i = 0; i < ret.W; i++) {
            Range r = get_base_range(BASE_COMP_B[seq[i]]);
            for (size_t j = 1; j < k; j++) {
                r = get_neighbor(r, BASE_COMP_B[seq[i+j]]);
            }
            for (auto f = r.first; f <= r.second; f++) {
                auto rid = mref_to_ref_id(fm_to_pac(f));
                auto g = rid_groups[rid];
                ret.set_genome(g, i);
            }
        }

        return ret;
    }

    using FwdRevCoords = std::pair< std::vector<i64>, std::vector<i64> >;

    //Returns all FM index coordinates which translate into reference 
    //coordinates that overlap the specified range
    FwdRevCoords range_to_fms(std::string ref_name, i64 start, i64 end) {

        std::vector<i64> fwd_fms, rev_fms;

        auto ref_len = static_cast<i64>(size() / 2);

        auto slop = static_cast<int>( ceil(log(ref_len) / log(4)) );

        auto pac_min = ref_to_pac(ref_name, start),
             pac_max = pac_min + (end - start) - 1;

        i64 fwd_st;
        if (ref_len - pac_max > slop) {
            fwd_st = pac_max + slop;
        } else {
            fwd_st = ref_len - 1;
        }

        Range r = get_base_range(get_base(fwd_st));
        for (auto i = fwd_st-1; i >= pac_max && i <= fwd_st; i--) {
            r = get_neighbor(r, get_base(i));
        }

        for (auto f = r.first; f <= r.second; f++) {
            auto loc = fm_to_pac(f);
            if (loc == pac_max) {
                r = Range(f,f);
                break;
            }
        }

        fwd_fms.push_back(r.first);
        for (auto i = pac_max-1; i >= pac_min && i < pac_max; i--) {
            r = get_neighbor(r, get_base(i));
            fwd_fms.push_back(r.first);
        }

        i64 rev_st;
        if (pac_min > slop) {
            rev_st = pac_min - slop;
        } else {
            rev_st = 0;
        }

        r = get_base_range(BASE_COMP_B[get_base(rev_st)]);
        for (i64 i = rev_st+1; i <= pac_min; i++) {
            r = get_neighbor(r, BASE_COMP_B[get_base(i)]);
        }

        for (auto f = r.first; f <= r.second; f++) {
            auto loc = fm_to_mref(f);
            if (loc == pac_min) {
                r = Range(f,f);
                break;
            }
        }

        rev_fms.push_back(r.first);
        for (auto i = pac_min+1; i <= pac_max; i++) {
            r = get_neighbor(r, BASE_COMP_B[get_base(i)]);
            rev_fms.push_back(r.first);
        }

        return FwdRevCoords(rev_fms, fwd_fms);
    }
    
    private:
    bwt_t *index_;
    bntseq_t *bns_;
    u8 *pacseq_;
    bool loaded_;
    size_t size_;
};


#endif
