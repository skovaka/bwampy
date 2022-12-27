#include "bwampy.hpp"

#define PY_BWA_INDEX_METH(P) c.def(#P, &RefIndex::P);
#define PY_BWA_INDEX_VEC(P) c.def(#P, py::vectorize(&RefIndex::P));
namespace py = pybind11;

PYBIND11_MODULE(bwampy, m) {
    m.doc() = R"pbdoc(Python wrapper for BWA low-level FM-index mapping functions)pbdoc";

    py::class_<PanKmerBitvec> (m, "PanKmerBitvec", py::buffer_protocol())
        .def_buffer([](PanKmerBitvec &r) -> py::buffer_info {
             return py::buffer_info(
                r.vec.data(),
                sizeof(u64),
                py::format_descriptor<u64>::format(),
                2,
                {r.H, r.W},
                {sizeof(u64) * r.W, sizeof(u64)}
            );
    });

    py::class_<RefIndex> c(m, "RefIndex");

    c.def(py::init<>());
    c.def(py::init<const std::string &>());
    c.def(py::init<const std::string &, bool>());
    c.def(py::init<const std::string &, bool, bool>());
    PY_BWA_INDEX_METH(pan_kmer_bitvec);
    PY_BWA_INDEX_METH(create);
    PY_BWA_INDEX_METH(load_index);
    PY_BWA_INDEX_METH(bwt_loaded);
    PY_BWA_INDEX_METH(load_pacseq);
    PY_BWA_INDEX_METH(destroy);
    PY_BWA_INDEX_METH(get_neighbor);
    PY_BWA_INDEX_METH(get_base_range);
    PY_BWA_INDEX_METH(fm_to_pac);
    PY_BWA_INDEX_METH(fm_to_mref);
    PY_BWA_INDEX_VEC(mref_to_pac);
    PY_BWA_INDEX_VEC(pac_to_ref);
    PY_BWA_INDEX_METH(ref_to_pac);
    PY_BWA_INDEX_METH(size);
    //PY_BWA_INDEX_METH(mrefs_to_ref_coord);
    PY_BWA_INDEX_METH(get_seqs);
    PY_BWA_INDEX_METH(pacseq_loaded);
    PY_BWA_INDEX_METH(get_sa_loc);
    PY_BWA_INDEX_METH(get_ref_name);
    PY_BWA_INDEX_METH(get_ref_len);
    PY_BWA_INDEX_METH(get_pac_shift);
    //PY_BWA_INDEX_METH(range_to_fms);
    PY_BWA_INDEX_METH(is_mref_fwd);
    PY_BWA_INDEX_METH(is_mref_flipped);
    PY_BWA_INDEX_VEC(get_base);
    c.def("pac_to_ref_id", static_cast<i32 (RefIndex::*)(i64)> (&RefIndex::pac_to_ref_id) );
    c.def("mref_to_ref_id", static_cast<i32 (RefIndex::*)(i64)> (&RefIndex::mref_to_ref_id) );
    c.def("get_ref_id", static_cast<i32 (RefIndex::*)(const std::string &)> (&RefIndex::get_ref_id));

}
