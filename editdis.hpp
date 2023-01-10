// Code Source: https://github.com/zhangxianbing/editdis
#pragma once

#include <algorithm>
#include <bitset>
#include <codecvt>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <locale>
#include <map>
#include <string>
#include <vector>

class editdis {
public:
    static inline int edit_distance(const std::string &a, const std::string &b,
                                    bool convert = true) {
        if (convert) {
            return edit_distance(to_wstring(a), to_wstring(b));
        }
        std::vector<int64_t> ia{a.begin(), a.end()}, ib{b.begin(), b.end()};
        return _edit_distance(ia.data(), ia.size(), ib.data(), ib.size());
    }

    static inline int edit_distance(const std::wstring &a, const std::wstring &b) {
        std::vector<int64_t> ia{a.begin(), a.end()}, ib{b.begin(), b.end()};
        return _edit_distance(ia.data(), ia.size(), ib.data(), ib.size());
    }

private:
    static std::wstring to_wstring(const std::string &s) {
        std::wstring_convert<std::codecvt_utf8<wchar_t>> converter;
        auto ws = converter.from_bytes(s);
        return ws;
    }

    // forking from:
    // https://github.com/roy-ht/editdistance/blob/master/editdistance/_editdistance.cpp
    template <typename T, typename TVALUE>
    static int _edit_distance_bpv(T &cmap, int64_t const *vec, size_t const &vecsize,
                                  unsigned int const &tmax, unsigned int const &tlen) {
        int D = tmax * 64 + tlen;
        TVALUE D0, HP, HN, VP, VN;
        uint64_t top = (1LL << (tlen - 1));
        uint64_t lmb = (1LL << 63);

        for (size_t i = 0; i <= tmax; ++i) {
            VP[i] = 0;
            VN[i] = 0;
        }
        for (size_t i = 0; i < tmax; ++i) VP[i] = ~0;
        for (size_t i = 0; i < tlen; ++i) VP[tmax] |= (1LL << i);
        for (size_t i = 0; i < vecsize; ++i) {
            TVALUE &PM = cmap[vec[i]];
            for (int r = 0; r <= tmax; ++r) {
                uint64_t X = PM[r];
                if (r > 0 && (HN[r - 1] & lmb)) X |= 1LL;
                D0[r] = (((X & VP[r]) + VP[r]) ^ VP[r]) | X | VN[r];
                HP[r] = VN[r] | ~(D0[r] | VP[r]);
                HN[r] = D0[r] & VP[r];
                X = (HP[r] << 1LL);
                if (r == 0 || HP[r - 1] & lmb) X |= 1LL;
                VP[r] = (HN[r] << 1LL) | ~(D0[r] | X);
                if (r > 0 && (HN[r - 1] & lmb)) VP[r] |= 1LL;
                VN[r] = D0[r] & X;
            }
            if (HP[tmax] & top)
                ++D;
            else if (HN[tmax] & top)
                --D;
        }
        return D;
    }

    /// c.f. http://handasse.blogspot.com/2009/04/c_29.html
    template <typename T>
    static int _edit_distance_dp(T const *str1, size_t const size1, T const *str2,
                                 size_t const size2) {
        std::vector<std::vector<uint32_t>> d(2, std::vector<uint32_t>(size2 + 1));
        d[0][0] = 0;
        d[1][0] = 1;
        for (int i = 0; i < size2 + 1; i++) d[0][i] = i;
        for (int i = 1; i < size1 + 1; i++) {
            for (int j = 1; j < size2 + 1; j++) {
                d[i & 1][j] =
                    std::min(std::min(d[(i - 1) & 1][j], d[i & 1][j - 1]) + 1,
                             d[(i - 1) & 1][j - 1] + (str1[i - 1] == str2[j - 1] ? 0 : 1));
            }
        }
        return d[size1 & 1][size2];
    }

    template <size_t N>
    struct varr {
        uint64_t arr_[N];
        uint64_t &operator[](size_t const &i) { return arr_[i]; }
    };

    template <size_t N>
    static int _edit_distance_map(int64_t const *a, size_t const asize, int64_t const *b,
                                  size_t const bsize) {
        typedef std::map<int64_t, varr<N>> cmap_v;
        cmap_v cmap;
        unsigned int tmax = (asize - 1) >> 6;
        unsigned int tlen = asize - tmax * 64;
        for (size_t i = 0; i < tmax; ++i) {
            for (size_t j = 0; j < 64; ++j) cmap[a[i * 64 + j]][i] |= (1LL << j);
        }
        for (size_t i = 0; i < tlen; ++i) cmap[a[tmax * 64 + i]][tmax] |= (1LL << i);
        return _edit_distance_bpv<cmap_v, typename cmap_v::mapped_type>(cmap, b, bsize, tmax,
                                                                        tlen);
    }

    static int _edit_distance(const int64_t *a, const unsigned int asize, const int64_t *b,
                              const unsigned int bsize) {
        if (asize == 0)
            return bsize;
        else if (bsize == 0)
            return asize;
        int64_t const *ap, *bp;
        unsigned int const *asizep, *bsizep;
        if (asize < bsize)
            ap = b, bp = a, asizep = &bsize, bsizep = &asize;
        else
            ap = a, bp = b, asizep = &asize, bsizep = &bsize;
        size_t vsize = ((*asizep - 1) >> 6) + 1;
        if (vsize > 10) {
            int64_t const *_ = ap;
            unsigned int const *__ = asizep;
            ap = bp, bp = _, asizep = bsizep, bsizep = __;
            vsize = ((*asizep - 1) >> 6) + 1;
        }

        if (vsize == 1)
            return _edit_distance_map<1>(ap, *asizep, bp, *bsizep);
        else if (vsize == 2)
            return _edit_distance_map<2>(ap, *asizep, bp, *bsizep);
        else if (vsize == 3)
            return _edit_distance_map<3>(ap, *asizep, bp, *bsizep);
        else if (vsize == 4)
            return _edit_distance_map<4>(ap, *asizep, bp, *bsizep);
        else if (vsize == 5)
            return _edit_distance_map<5>(ap, *asizep, bp, *bsizep);
        else if (vsize == 6)
            return _edit_distance_map<6>(ap, *asizep, bp, *bsizep);
        else if (vsize == 7)
            return _edit_distance_map<7>(ap, *asizep, bp, *bsizep);
        else if (vsize == 8)
            return _edit_distance_map<8>(ap, *asizep, bp, *bsizep);
        else if (vsize == 9)
            return _edit_distance_map<9>(ap, *asizep, bp, *bsizep);
        else if (vsize == 10)
            return _edit_distance_map<10>(ap, *asizep, bp, *bsizep);
        return _edit_distance_dp<int64_t>(ap, *asizep, bp, *bsizep);
    }
};
