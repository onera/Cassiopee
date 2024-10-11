#include "icapsule.h"
#include "common/Karray.h"
#include "point.h"
#include "io.h"

Smesh IMesh::make_patch(const E_Float *ptag)
{
    patch.clear();

    for (E_Int fid : skin) {
        const auto &pn = F[fid];
        bool is_patch = true;
        for (E_Int pid : pn) {
            if (ptag[pid] != 1) {
                is_patch = false;
                break;
            }
        }
        if (is_patch)
            patch.insert(fid);
    }
    
    return Smesh(*this);
}

Smesh IMesh::make_patch(const Smesh &spatch)
{
    //patch.clear();

    return Smesh();
}

ICapsule::ICapsule(const Karray &marray, const std::vector<Karray> &sarrays,
    const std::vector<E_Float *> &ptags)
{
    M = IMesh(marray);
    M.make_skin();
    M.orient_skin(OUT);
    M.triangulate_skin();
    M.make_bbox();
    M.hash_skin();

    Ss.reserve(sarrays.size());
    spatches.reserve(sarrays.size());

    std::vector<std::vector<PointLoc>> plocs;
    plocs.reserve(sarrays.size());

    for (size_t i = 0; i < sarrays.size(); i++) {
        Ss.push_back(IMesh(sarrays[i]));
        Ss[i].make_skin();
        Ss[i].orient_skin(IN);
        Ss[i].triangulate_skin();
        spatches.push_back(Ss[i].make_patch(ptags[i]));
        plocs.push_back(M.locate(spatches[i]));
    }

    // Correct near-vertex/edge situations
    std::vector<Point> points;
    for (size_t i = 0; i < sarrays.size(); i++) {
        const Smesh &Sf = spatches[i];
        const std::vector<PointLoc> &locs = plocs[i];
        assert(locs.size() == (size_t)spatches[i].np);
        E_Int on_vertex = 0, on_edge = 0;
        for (size_t j = 0; j < locs.size(); j++) {
            if (locs[j].v_idx != -1) on_vertex++;
            else if (locs[j].e_idx != -1) {
                on_edge++;
                points.push_back(Point(Sf.X[j], Sf.Y[j], Sf.Z[j]));
            }
        }
        printf("on vertex: %d - on edge: %d\n", on_vertex, on_edge);
    }

    point_write("points", points); 
    
    mpatches.reserve(sarrays.size());
}


PyObject *K_XCORE::icapsule_extract_master(PyObject *self, PyObject *args)
{
    PyObject *ICAPSULE;
    if (!PYPARSETUPLE_(args, O_, &ICAPSULE)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(ICAPSULE, "ICapsule")) {
        RAISE("Bad ICapsule hook.");
        return NULL;
    }

    ICapsule *icap = (ICapsule *)PyCapsule_GetPointer(ICAPSULE, "ICapsule");

    auto Mout = icap->M.export_karray();

    return Mout;
}

PyObject *K_XCORE::icapsule_extract_slave(PyObject *self, PyObject *args)
{
    PyObject *ICAPSULE;
    E_Int INDEX;
    if (!PYPARSETUPLE_(args, O_ I_, &ICAPSULE, &INDEX)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(ICAPSULE, "ICapsule")) {
        RAISE("Bad ICapsule hook.");
        return NULL;
    }

    ICapsule *icap = (ICapsule *)PyCapsule_GetPointer(ICAPSULE, "ICapsule");

    auto Sout = icap->Ss[INDEX].export_karray();

    return Sout;
}

PyObject *K_XCORE::icapsule_init(PyObject *self, PyObject *args)
{
    PyObject *MASTER, *SLAVES, *PTAGS;
    if (!PYPARSETUPLE_(args, OOO_, &MASTER, &SLAVES, &PTAGS)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray marray;
    if (Karray_parse_ngon(MASTER, marray) != 0) {
        RAISE("Master should be an NGon.");
        return NULL;
    }

    E_Int nslaves = PyList_Size(SLAVES);
    std::vector<Karray> sarrays(nslaves); 
    E_Int ok = 0;

    for (E_Int i = 0; i < nslaves; i++) {
        PyObject *SLAVE = PyList_GetItem(SLAVES, i);
        E_Int ret = Karray_parse_ngon(SLAVE, sarrays[i]);
        Py_DECREF(SLAVE);
        if (ret != 0) {
            RAISE("Slaves should be NGons.");
            break;
        }
        ok++;
    }

    if (ok != nslaves) {
        Karray_free_ngon(marray);
        for (E_Int i = 0; i < ok; i++)
            Karray_free_ngon(sarrays[i]);
        return NULL;
    }

    std::vector<E_Float *> ptags(nslaves, NULL);
    
    if (PyList_Size(PTAGS) != nslaves) {
        RAISE("Ptags should be the same size as slaves.");
        Karray_free_ngon(marray);
        for (E_Int i = 0; i < nslaves; i++)
            Karray_free_ngon(sarrays[i]);
        return NULL;
    }

    for (E_Int i = 0; i < nslaves; i++) {
        PyObject *PTAG = PyList_GetItem(PTAGS, i);
        E_Int size = -1;
        E_Int ret = K_NUMPY::getFromNumpyArray(PTAG, ptags[i], size, true);
        Py_DECREF(PTAG);
        if (ret != 1 || size != sarrays[i].npoints()) {
            RAISE("Ptag[i] should have size sarrays[i].npoints.");
            Karray_free_ngon(marray);
            for (E_Int j = 0; j < nslaves; j++)
                Karray_free_ngon(sarrays[j]);
            return NULL;
        }
    }

    ICapsule *icap = new ICapsule(marray, sarrays, ptags);

    PyObject *out = PyCapsule_New((void *)icap, "ICapsule", NULL);

    Karray_free_ngon(marray);
    for (E_Int i = 0; i < nslaves; i++)
        Karray_free_ngon(sarrays[i]);

    return out;
}
