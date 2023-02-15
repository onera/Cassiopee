/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr); Imad Hammani (imad.hammani@onera.fr)

#ifndef NUGA_METRIC_SENSOR_HXX
#define NUGA_METRIC_SENSOR_HXX

#include "Nuga/include/sensor.hxx"
#include "Nuga/include/Metric.h"
//#include "Nuga/include/shell_smoother.hxx"
//#include "Nuga/include/V1_smoother.hxx"


namespace NUGA
{

template <typename mesh_t>
class metric_sensor : public sensor<mesh_t, DELAUNAY::VarMetric<DELAUNAY::Aniso3D>>
{
  public:
    bool _single_pass_done;
    
    //using cell_incr_t = DELAUNAY::Aniso3D;
    using sensor_input_t = DELAUNAY::VarMetric<DELAUNAY::Aniso3D>; // either vector<int> (ISO) or vector(incr_t<3>) (DIR)
    using parent_t = sensor<mesh_t, sensor_input_t>;
    using output_t = typename mesh_t::output_t; //fixme: static assert to add : must be ISO => IntVec

    metric_sensor(mesh_t& mesh) : parent_t(mesh, nullptr)
    {
        parent_t::_hmesh._sensor = this; /*in order to have an access to the field values from hmesh*/
        _canon_info.resize(mesh._ng.PHs.size());
        _start_nodes.resize(mesh._ng.PHs.size());
    }

    E_Int assign_data(const K_FLD::FloatArray& data);
    
    bool fill_adap_incr(output_t& adap_incr, bool do_agglo) override;

    void metric_fix(mesh_t& hmesh, output_t& adap_incr);

    bool stop() override { return _single_pass_done; }

    void Q4_adap_compute(E_Int PGi, output_t& adap_incr);

    void Hexa_adap_compute(E_Int PHi, output_t& adap_incr);

    bool is_FXY(E_Int PGi, output_t& adap_incr);

    bool is_HX12(E_Int PHi, output_t& adap_incr);

    bool is_HX18(E_Int PHi, output_t& adap_incr);

    bool is_HX27(E_Int PHi, output_t& adap_incr);

    std::vector<std::vector<E_Int>> _canon_info;
    std::vector<std::vector<E_Int>> _start_nodes;

    void compute_canon_info_bottom(E_Int PHi, E_Int n0);
    void compute_canon_info_top(E_Int PHi, E_Int n0);
    void compute_canon_info_left(E_Int PHi, E_Int n0);
    void compute_canon_info_right(E_Int PHi, E_Int n0);
    void compute_canon_info_front(E_Int PHi, E_Int n0);
    void compute_canon_info_back(E_Int PHi, E_Int n0);

    bool FIX(E_Int PHi, E_Int nei, output_t& adap_incr);
};

template <typename mesh_t>
void metric_sensor<mesh_t>::compute_canon_info_bottom(E_Int PHi, E_Int n0)
{
    auto& hmesh = parent_t::_hmesh;
    const auto& ng = hmesh._ng;
    const auto& F2E = hmesh._F2E;

    E_Int nodes[8];

    const E_Int *faces = ng.PHs.get_facets_ptr(PHi);

    bool reorient;
    E_Int i0;
    E_Int local[4];
    E_Int canon_data[12];
    E_Int PGi;

    _start_nodes[PHi].resize(6);

    // Bottom
    PGi = faces[0] - 1;
    const E_Int *pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = n0;
    canon_data[0] = (E_Int) reorient;
    canon_data[1] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    nodes[0] = local[0];
    nodes[1] = local[1];
    nodes[2] = local[2];
    nodes[3] = local[3];
    _start_nodes[PHi][0] = i0;

    // Left
    PGi = faces[2] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[1] == nodes[3]);
    nodes[7] = local[2];
    nodes[4] = local[3];
    canon_data[4] = (E_Int) reorient;
    canon_data[5] = i0;
    _start_nodes[PHi][2] = i0;

    // Right
    PGi = faces[3] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[2]);
    nodes[6] = local[2];
    nodes[5] = local[3];
    canon_data[6] = (E_Int) reorient;
    canon_data[7] = i0;
    _start_nodes[PHi][3] = i0;

    // Top
    PGi = faces[1] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[4]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[4]);
    assert(local[1] == nodes[5]);
    assert(local[2] == nodes[6]);
    assert(local[3] == nodes[7]);
    canon_data[2] = (E_Int) reorient;
    canon_data[3] = i0;
    _start_nodes[PHi][1] = i0;

    // Front
    PGi = faces[4] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[0]);
    assert(local[2] == nodes[4]);
    assert(local[3] == nodes[5]);
    canon_data[8] = (E_Int) reorient;
    canon_data[9] = i0;
    _start_nodes[PHi][4] = i0;

    // Back
    PGi = faces[5] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[2]);
    assert(local[1] == nodes[3]);
    assert(local[2] == nodes[7]);
    assert(local[3] == nodes[6]);
    canon_data[10] = (E_Int) reorient;
    canon_data[11] = i0;
    _start_nodes[PHi][5] = i0;

    _canon_info[PHi].resize(6);

    // deduce swap
    for (E_Int i = 0; i < 6; i++) {
        E_Int reorient = canon_data[2*i];
        E_Int i0 = canon_data[2*i+1];
        if (reorient) {
            if (i0%2 == 0) _canon_info[PHi][i] = 1;
            else _canon_info[PHi][i] = 0;
        } else {
          if (i0%2 == 0) _canon_info[PHi][i] = 0;
          else _canon_info[PHi][i] = 1;
        }
    }
}

template <typename mesh_t>
void metric_sensor<mesh_t>::compute_canon_info_top(E_Int PHi, E_Int n0)
{
    auto& hmesh = parent_t::_hmesh;
    const auto& ng = hmesh._ng;
    const auto& F2E = hmesh._F2E;

    E_Int nodes[8];

    const E_Int *faces = ng.PHs.get_facets_ptr(PHi);

    bool reorient;
    E_Int i0;
    E_Int local[4];
    E_Int canon_data[12];
    E_Int PGi;

    _start_nodes[PHi].resize(6);

    // Top
    PGi = faces[1] - 1;
    const E_Int *pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = n0;
    canon_data[2] = (E_Int) reorient;
    canon_data[3] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == pN[i0]);
    nodes[4] = local[0];
    nodes[5] = local[1];
    nodes[6] = local[2];
    nodes[7] = local[3];
    _start_nodes[PHi][1] = i0;

    // Left
    PGi = faces[2] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[7]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[7]);
    assert(local[1] == nodes[4]);
    nodes[0] = local[2];
    nodes[3] = local[3];
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, local[2]);
    canon_data[4] = (E_Int) reorient;
    canon_data[5] = i0;
    _start_nodes[PHi][2] = i0;

    // Right
    PGi = faces[3] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[6]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[6]);
    assert(local[1] == nodes[5]);
    nodes[1] = local[2];
    nodes[2] = local[3];
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, local[2]);
    canon_data[6] = (E_Int) reorient;
    canon_data[7] = i0;
    _start_nodes[PHi][3] = i0;

    // bottom
    PGi = faces[0] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[1] == nodes[1]);
    assert(local[2] == nodes[2]);
    assert(local[3] == nodes[3]);
    canon_data[0] = (E_Int) reorient;
    canon_data[1] = i0;
    _start_nodes[PHi][0] = i0;

    // Front
    PGi = faces[4] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    canon_data[8] = (E_Int) reorient;
    canon_data[9] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[0]);
    assert(local[2] == nodes[4]);
    assert(local[3] == nodes[5]);
    _start_nodes[PHi][4] = i0;

    // Back
    PGi = faces[5] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
    canon_data[10] = (E_Int) reorient;
    canon_data[11] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[2]);
    assert(local[1] == nodes[3]);
    assert(local[2] == nodes[7]);
    assert(local[3] == nodes[6]);
    _start_nodes[PHi][5] = i0;

    _canon_info[PHi].resize(6);

    // deduce swap
    for (E_Int i = 0; i < 6; i++) {
        E_Int reorient = canon_data[2*i];
        E_Int i0 = canon_data[2*i+1];
        if (reorient) {
            if (i0%2 == 0) _canon_info[PHi][i] = 1;
            else _canon_info[PHi][i] = 0;
        } else {
          if (i0%2 == 0) _canon_info[PHi][i] = 0;
          else _canon_info[PHi][i] = 1;
        }
    }
}

template <typename mesh_t>
void metric_sensor<mesh_t>::compute_canon_info_left(E_Int PHi, E_Int n0)
{
    auto& hmesh = parent_t::_hmesh;
    const auto& ng = hmesh._ng;
    const auto& F2E = hmesh._F2E;

    E_Int nodes[8];

    const E_Int *faces = ng.PHs.get_facets_ptr(PHi);

    bool reorient;
    E_Int i0;
    E_Int local[4];
    E_Int canon_data[12];
    E_Int PGi;

    _start_nodes[PHi].resize(6);

    // Left
    PGi = faces[2] - 1;
    const E_Int *pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = n0;
    canon_data[4] = (E_Int) reorient;
    canon_data[5] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == pN[i0]);
    nodes[0] = local[0];
    nodes[3] = local[1];
    nodes[7] = local[2];
    nodes[4] = local[3];
    _start_nodes[PHi][2] = i0;

    // Bottom
    PGi = faces[0] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    canon_data[0] = (E_Int) reorient;
    canon_data[1] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[3] == nodes[3]);
    nodes[1] = local[1];
    nodes[2] = local[2];
    _start_nodes[PHi][0] = i0;

    // Top
    PGi = faces[1] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[4]);
    canon_data[2] = (E_Int) reorient;
    canon_data[3] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[4]);
    assert(local[3] == nodes[7]);
    nodes[5] = local[1];
    nodes[6] = local[2];
    _start_nodes[PHi][1] = i0;

    // Right
    PGi = faces[3] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[2]);
    assert(local[2] == nodes[6]);
    assert(local[3] == nodes[5]);
    canon_data[6] = (E_Int) reorient;
    canon_data[7] = i0;
    _start_nodes[PHi][3] = i0;

    // Front
    PGi = faces[4] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    canon_data[8] = (E_Int) reorient;
    canon_data[9] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[0]);
    assert(local[2] == nodes[4]);
    assert(local[3] == nodes[5]);
    _start_nodes[PHi][4] = i0;

    // Back
    PGi = faces[5] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
    canon_data[10] = (E_Int) reorient;
    canon_data[11] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[2]);
    assert(local[1] == nodes[3]);
    assert(local[2] == nodes[7]);
    assert(local[3] == nodes[6]);
    _start_nodes[PHi][5] = i0;

    _canon_info[PHi].resize(6);

    // deduce swap
    for (E_Int i = 0; i < 6; i++) {
        E_Int reorient = canon_data[2*i];
        E_Int i0 = canon_data[2*i+1];
        if (reorient) {
            if (i0%2 == 0) _canon_info[PHi][i] = 1;
            else _canon_info[PHi][i] = 0;
        } else {
          if (i0%2 == 0) _canon_info[PHi][i] = 0;
          else _canon_info[PHi][i] = 1;
        }
    }
}

template <typename mesh_t>
void metric_sensor<mesh_t>::compute_canon_info_right(E_Int PHi, E_Int n0)
{
    auto& hmesh = parent_t::_hmesh;
    const auto& ng = hmesh._ng;
    const auto& F2E = hmesh._F2E;

    E_Int nodes[8];

    const E_Int *faces = ng.PHs.get_facets_ptr(PHi);

    bool reorient;
    E_Int i0;
    E_Int local[4];
    E_Int canon_data[12];
    E_Int PGi;

    _start_nodes[PHi].resize(6);

    // Right
    PGi = faces[3] - 1;
    const E_Int *pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = n0;
    canon_data[6] = (E_Int) reorient;
    canon_data[7] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == pN[i0]);
    nodes[1] = local[0];
    nodes[2] = local[1];
    nodes[6] = local[2];
    nodes[5] = local[3];
    _start_nodes[PHi][3] = i0;

    // Front
    PGi = faces[4] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    canon_data[8] = (E_Int) reorient;
    canon_data[9] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[3] == nodes[5]);
    nodes[0] = local[1];
    nodes[4] = local[2];
    _start_nodes[PHi][4] = i0;

    // Back
    PGi = faces[5] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
    canon_data[10] = (E_Int) reorient;
    canon_data[11] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[2]);
    assert(local[3] == nodes[6]);
    nodes[3] = local[1];
    nodes[7] = local[2];
    _start_nodes[PHi][5] = i0;

    // Left
    PGi = faces[2] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[1] == nodes[3]);
    assert(local[2] == nodes[7]);
    assert(local[3] == nodes[4]);
    canon_data[4] = (E_Int) reorient;
    canon_data[5] = i0;
    _start_nodes[PHi][2] = i0;

    // Bottom
    PGi = faces[0] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    canon_data[0] = (E_Int) reorient;
    canon_data[1] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[1] == nodes[1]);
    assert(local[2] == nodes[2]);
    assert(local[3] == nodes[3]);
    _start_nodes[PHi][0] = i0;

    // Top
    PGi = faces[1] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[4]);
    canon_data[2] = (E_Int) reorient;
    canon_data[3] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[4]);
    assert(local[1] == nodes[5]);
    assert(local[2] == nodes[6]);
    assert(local[3] == nodes[7]);
    _start_nodes[PHi][1] = i0;

    _canon_info[PHi].resize(6);

    // deduce swap
    for (E_Int i = 0; i < 6; i++) {
        E_Int reorient = canon_data[2*i];
        E_Int i0 = canon_data[2*i+1];
        if (reorient) {
            if (i0%2 == 0) _canon_info[PHi][i] = 1;
            else _canon_info[PHi][i] = 0;
        } else {
          if (i0%2 == 0) _canon_info[PHi][i] = 0;
          else _canon_info[PHi][i] = 1;
        }
    }
}

template <typename mesh_t>
void metric_sensor<mesh_t>::compute_canon_info_front(E_Int PHi, E_Int n0)
{
    auto& hmesh = parent_t::_hmesh;
    const auto& ng = hmesh._ng;
    const auto& F2E = hmesh._F2E;

    E_Int nodes[8];

    const E_Int *faces = ng.PHs.get_facets_ptr(PHi);

    bool reorient;
    E_Int i0;
    E_Int local[4];
    E_Int canon_data[12];
    E_Int PGi;

    _start_nodes[PHi].resize(6);

    // Front
    PGi = faces[4] - 1;
    const E_Int *pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = n0;
    canon_data[8] = (E_Int) reorient;
    canon_data[9] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == pN[i0]);
    nodes[1] = local[0];
    nodes[0] = local[1];
    nodes[4] = local[2];
    nodes[5] = local[3];
    _start_nodes[PHi][4] = i0;

    // Bottom
    PGi = faces[0] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    canon_data[0] = (E_Int) reorient;
    canon_data[1] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[1] == nodes[1]);
    nodes[2] = local[2];
    nodes[3] = local[3];
    _start_nodes[PHi][0] = i0;

    // Top
    PGi = faces[1] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[4]);
    canon_data[2] = (E_Int) reorient;
    canon_data[3] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[4]);
    assert(local[1] == nodes[5]);
    nodes[6] = local[2];
    nodes[7] = local[3];
    _start_nodes[PHi][1] = i0;

    // Right
    PGi = faces[3] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[2]);
    assert(local[2] == nodes[6]);
    assert(local[3] == nodes[5]);
    canon_data[6] = (E_Int) reorient;
    canon_data[7] = i0;
    _start_nodes[PHi][3] = i0;

    // Left
    PGi = faces[2] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    canon_data[4] = (E_Int) reorient;
    canon_data[5] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[1] == nodes[3]);
    assert(local[2] == nodes[7]);
    assert(local[3] == nodes[4]);
    _start_nodes[PHi][2] = i0;

    // Back
    PGi = faces[5] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
    canon_data[10] = (E_Int) reorient;
    canon_data[11] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[2]);
    assert(local[1] == nodes[3]);
    assert(local[2] == nodes[7]);
    assert(local[3] == nodes[6]);
    _start_nodes[PHi][5] = i0;

    _canon_info[PHi].resize(6);

    // deduce swap
    for (E_Int i = 0; i < 6; i++) {
        E_Int reorient = canon_data[2*i];
        E_Int i0 = canon_data[2*i+1];
        if (reorient) {
            if (i0%2 == 0) _canon_info[PHi][i] = 1;
            else _canon_info[PHi][i] = 0;
        } else {
          if (i0%2 == 0) _canon_info[PHi][i] = 0;
          else _canon_info[PHi][i] = 1;
        }
    }
}

template <typename mesh_t>
void metric_sensor<mesh_t>::compute_canon_info_back(E_Int PHi, E_Int n0)
{
    auto& hmesh = parent_t::_hmesh;
    const auto& ng = hmesh._ng;
    const auto& F2E = hmesh._F2E;

    E_Int nodes[8];

    const E_Int *faces = ng.PHs.get_facets_ptr(PHi);

    bool reorient;
    E_Int i0;
    E_Int local[4];
    E_Int canon_data[12];
    E_Int PGi;

    _start_nodes[PHi].resize(6);

    // Back
    PGi = faces[5] - 1;
    const E_Int *pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = n0;
    canon_data[10] = (E_Int) reorient;
    canon_data[11] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == pN[i0]);
    nodes[2] = local[0];
    nodes[3] = local[1];
    nodes[7] = local[2];
    nodes[6] = local[3];
    _start_nodes[PHi][5] = i0;

    // Bottom
    PGi = faces[0] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[2]);
    assert(local[1] == nodes[3]);
    nodes[0] = local[2];
    nodes[1] = local[3];
    // local[2] is the usual i0 here
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, local[2]);
    canon_data[0] = (E_Int) reorient;
    canon_data[1] = i0;
    _start_nodes[PHi][0] = i0;

    // Top
    PGi = faces[1] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[6]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[6]);
    assert(local[1] == nodes[7]);
    nodes[4] = local[2];
    nodes[5] = local[3];
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, local[2]);
    canon_data[2] = (E_Int) reorient;
    canon_data[3] = i0;
    _start_nodes[PHi][1] = i0;

    // Right
    PGi = faces[3] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[2]);
    assert(local[2] == nodes[6]);
    assert(local[3] == nodes[5]);
    canon_data[6] = (E_Int) reorient;
    canon_data[7] = i0;
    _start_nodes[PHi][3] = i0;

    // Left
    PGi = faces[2] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    canon_data[4] = (E_Int) reorient;
    canon_data[5] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[1] == nodes[3]);
    assert(local[2] == nodes[7]);
    assert(local[3] == nodes[4]);
    _start_nodes[PHi][2] = i0;

    // Front
    PGi = faces[4] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    canon_data[8] = (E_Int) reorient;
    canon_data[9] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[0]);
    assert(local[2] == nodes[4]);
    assert(local[3] == nodes[5]);
    _start_nodes[PHi][4] = i0;

    _canon_info[PHi].resize(6);

    // deduce swap
    for (E_Int i = 0; i < 6; i++) {
        E_Int reorient = canon_data[2*i];
        E_Int i0 = canon_data[2*i+1];
        if (reorient) {
            if (i0%2 == 0) _canon_info[PHi][i] = 1;
            else _canon_info[PHi][i] = 0;
        } else {
          if (i0%2 == 0) _canon_info[PHi][i] = 0;
          else _canon_info[PHi][i] = 1;
        }
    }
}

template <typename mesh_t>
E_Int metric_sensor<mesh_t>::assign_data(const K_FLD::FloatArray& arr) {
    // convert data to VarMetric<Aniso3D>

    E_Int npts = parent_t::_hmesh._crd.cols();

    assert(npts == arr.cols()); // input must be sized as number of pts (if decided to do differently,  need here to resize and fill missing field values)

    sensor_input_t data(parent_t::_hmesh._crd, -1, -1);
    DELAUNAY::Aniso3D m;

    // populate _data
    for (E_Int i = 0; i < npts; i++) {
        for (E_Int j = 0; j < 6; j++) {
            m[j] = arr(j,i);
        }
        if (!data.isValidMetric(m)) {
            fprintf(stderr, "BAD METRIC!\n");
            std::cout << "Point: " << i << std::endl;
            std::cout << m << std::endl;
            exit(1);
        }
        //assert(data.isValidMetric(m));
        data.setMetric(i,m);
    }

#ifdef DEBUG_METRIC_SENSOR
    std::cout << "assign_data ok\n";
#endif

    parent_t::_data = data;

  _single_pass_done = false;//reinit to enable

  return 0;

}


template <typename mesh_t>
bool metric_sensor<mesh_t>::is_FXY(E_Int PGi, output_t& adap_incr)
{
    if (adap_incr.face_adap_incr[PGi].n[0] > 0 && adap_incr.face_adap_incr[PGi].n[1] > 0)
        return true;
    return false;
}

template <typename mesh_t>
bool metric_sensor<mesh_t>::is_HX12(E_Int PHi, output_t& adap_incr)
{
    E_Int count = 0;
    for (E_Int i = 0; i < 3; i++) {
        if (adap_incr.cell_adap_incr[PHi].n[i] > 0) count++;
    }

    return (count == 1);
}

template <typename mesh_t>
bool metric_sensor<mesh_t>::is_HX18(E_Int PHi, output_t& adap_incr)
{
    E_Int count = 0;
    for (E_Int i = 0; i < 3; i++) {
        if (adap_incr.cell_adap_incr[PHi].n[i] > 0) count++;
    }

    return (count == 2);
}

template <typename mesh_t>
bool metric_sensor<mesh_t>::is_HX27(E_Int PHi, output_t& adap_incr)
{
    E_Int count = 0;
    for (E_Int i = 0; i < 3; i++) {
        if (adap_incr.cell_adap_incr[PHi].n[i] > 0) count++;
    }

    return (count == 3);
}

template <typename mesh_t>
bool metric_sensor<mesh_t>::FIX(E_Int PHi, E_Int nei, output_t& adap_incr)
{
    if (nei == IDX_NONE) return false;

    // both cells are in coherent config
    E_Int ao = adap_incr.cell_adap_incr[PHi].n[0];
    E_Int bo = adap_incr.cell_adap_incr[PHi].n[1];
    E_Int co = adap_incr.cell_adap_incr[PHi].n[2];
    E_Int an = adap_incr.cell_adap_incr[nei].n[0];
    E_Int bn = adap_incr.cell_adap_incr[nei].n[1];
    E_Int cn = adap_incr.cell_adap_incr[nei].n[2];

    E_Int PHi_to_mod = -1;

    // fix alphas
    if (::abs(ao-an) <= 1) goto fix_beta;
    PHi_to_mod = ao > an ? nei : PHi;
    adap_incr.cell_adap_incr[PHi_to_mod].n[0] = std::max(ao, an) - 1;

    // fix betas
    fix_beta:
    if (::abs(bo-bn) <= 1) goto fix_gamma;
    PHi_to_mod = bo > bn ? nei : PHi;
    adap_incr.cell_adap_incr[PHi_to_mod].n[1] = std::max(bo, bn) - 1;

    // fix gammas
    fix_gamma:
    if (::abs(co-cn) <= 1) goto end;
    PHi_to_mod = co > cn ? nei : PHi;
    adap_incr.cell_adap_incr[PHi_to_mod].n[2] = std::max(co, cn) - 1;
    
    end:
    bool has_changed = !(adap_incr.cell_adap_incr[PHi].n[0] == ao &&
             adap_incr.cell_adap_incr[PHi].n[1] == bo &&
             adap_incr.cell_adap_incr[PHi].n[2] == co);

    return has_changed;
}


/*
template <typename mesh_t>
void metric_sensor<mesh_t>::fix_adap_incr(mesh_t& hmesh, output_t& adap_incr)
{
    using cell_incr_t = typename output_t::cell_incr_t;
    using face_incr_t = typename output_t::face_incr_t;

    std::stack<E_Int> stk;

    for (E_Int PGi = 0; PGi < hmesh._ng.PGs.size(); PGi++) {
        if (adap_incr.face_adap_incr[PGi] != 0) stk.push(PGi);
    }

    while (!stk.empty()) {
        E_Int PGi = stk.top();
        stk.pop();

        E_Int own = hmesh._F2E(0, PGi);
        E_Int nei = hmesh._F2E(1, PGi);
        bool ok_nei = (nei != IDX_NONE && hmesh._PHtree.is_enabled(nei));
        if (!ok_nei) continue;

        // PGi is right of own and left of nei
        E_Int *pOwn = hmesh._ng.PHs.get_facets_ptr(own);
        E_Int *pNei = hmesh._ng.PHs.get_facets_ptr(nei);

        while (pOwn[3] != PGi+1) K_CONNECT::IdTool::right_shift<6>(pOwn, 1);
        while (pNei[2] != PGi+1) K_CONNECT::IdTool::right_shift<6>(pNei, 1);
        
        K_MESH::Hexahedron::reorder_pgs_right(hmesh._ng, hmesh._F2E, own, 0);
        compute_canon_info_right(own, 0);
        E_Int i0 = _start_nodes[own][3];
        Hexa_adap_compute(own, adap_incr);

        K_MESH::Hexahedron::reorder_pgs_left(hmesh._ng, hmesh._F2E, nei, i0);
        compute_canon_info_left(nei, i0);
        Hexa_adap_compute(nei, adap_incr);

        E_Int ao = adap_incr.cell_adap_incr[own].n[0];
        E_Int bo = adap_incr.cell_adap_incr[own].n[1];
        E_Int co = adap_incr.cell_adap_incr[own].n[2];
        E_Int an = adap_incr.cell_adap_incr[nei].n[0];
        E_Int bn = adap_incr.cell_adap_incr[nei].n[1];
        E_Int cn = adap_incr.cell_adap_incr[nei].n[2];

        E_Int PHi_to_mod = -1;
        bool own_mod_a = false;
        bool nei_mod_a = false;
        bool own_mod_b = false;
        bool nei_mod_b = false;
        bool own_mod_c = false;
        bool nei_mod_c = false;

        // fix alphas
        if (::abs(ao-an) <= 1) goto fix_beta;
        PHi_to_mod = ao > an ? nei : own;
        if (PHi_to_mod == own) own_mod_a = true;
        if (PHi_to_mod == nei) nei_mod_a = true;
        adap_incr.cell_adap_incr[PHi_to_mod].n[0] = std::max(ao, an) - 1;

        // fix betas
        fix_beta:
        if (::abs(bo-bn) <= 1) goto fix_gamma;
        PHi_to_mod = bo > bn ? nei : own;
        if (PHi_to_mod == own) own_mod_b = true;
        if (PHi_to_mod == nei) nei_mod_b = true;
        adap_incr.cell_adap_incr[PHi_to_mod].n[1] = std::max(bo, bn) - 1;

        // fix gammas
        fix_gamma:
        if (::abs(co-cn) <= 1) goto update_faces;
        PHi_to_mod = co > cn ? nei : own;
        if (PHi_to_mod == own) own_mod_c = true;
        if (PHi_to_mod == nei) nei_mod_c = true;
        adap_incr.cell_adap_incr[PHi_to_mod].n[2] = std::max(co, cn) - 1;

        update_faces:
        // update faces and push to stack accordingly
        E_Int BOT = pOwn[0]-1;
        E_Int TOP = pOwn[1]-1;
        E_Int LFT = pOwn[2]-1;
        E_Int RGT = pOwn[3]-1;
        E_Int FRO = pOwn[4]-1;
        E_Int BCK = pOwn[5]-1;
        auto& bot_inc_own = adap_incr.face_adap_incr[BOT];
        auto& top_inc_own = adap_incr.face_adap_incr[TOP];
        auto& lft_inc_own = adap_incr.face_adap_incr[LFT];
        auto& rgt_inc_own = adap_incr.face_adap_incr[RGT];
        auto& fro_inc_own = adap_incr.face_adap_incr[FRO];
        auto& bck_inc_own = adap_incr.face_adap_incr[BCK];

        const auto& swap_own = _canon_info[own];

        ao = adap_incr.cell_adap_incr[own].n[0];
        bo = adap_incr.cell_adap_incr[own].n[1];
        co = adap_incr.cell_adap_incr[own].n[2];
        
        if (own_mod_a) {
            if (swap_own[0]) bot_inc_own.n[1] = std::max(bot_inc_own.n[1], ao);
            else             bot_inc_own.n[0] = std::max(bot_inc_own.n[0], ao);
            if (swap_own[1]) top_inc_own.n[1] = std::max(top_inc_own.n[1], ao);
            else             top_inc_own.n[0] = std::max(top_inc_own.n[0], ao);
            if (swap_own[4]) fro_inc_own.n[1] = std::max(fro_inc_own.n[1], ao);
            else             fro_inc_own.n[0] = std::max(fro_inc_own.n[0], ao);
            if (swap_own[5]) bck_inc_own.n[1] = std::max(bck_inc_own.n[1], ao);
            else             bck_inc_own.n[0] = std::max(bck_inc_own.n[0], ao);
            stk.push(BOT);
            stk.push(TOP);
            stk.push(FRO);
            stk.push(BCK);
        }

        if (own_mod_b) {
            if (swap_own[0]) bot_inc_own.n[0] = std::max(bot_inc_own.n[0], bo);
            else             bot_inc_own.n[1] = std::max(bot_inc_own.n[1], bo);
            if (swap_own[1]) top_inc_own.n[0] = std::max(top_inc_own.n[0], bo);
            else             top_inc_own.n[1] = std::max(top_inc_own.n[1], bo);
            if (swap_own[2]) lft_inc_own.n[1] = std::max(lft_inc_own.n[1], bo);
            else             lft_inc_own.n[0] = std::max(lft_inc_own.n[0], bo);
            if (swap_own[3]) rgt_inc_own.n[1] = std::max(rgt_inc_own.n[1], bo);
            else             rgt_inc_own.n[0] = std::max(rgt_inc_own.n[0], bo);
            stk.push(BOT);
            stk.push(TOP);
            stk.push(LFT);
        }

        if (own_mod_c) {
            if (swap_own[2]) lft_inc_own.n[0] = std::max(lft_inc_own.n[0], co);
            else             lft_inc_own.n[1] = std::max(lft_inc_own.n[1], co);
            if (swap_own[3]) rgt_inc_own.n[0] = std::max(rgt_inc_own.n[0], co);
            else             rgt_inc_own.n[1] = std::max(rgt_inc_own.n[1], co);
            if (swap_own[4]) fro_inc_own.n[0] = std::max(fro_inc_own.n[0], co);
            else             fro_inc_own.n[1] = std::max(fro_inc_own.n[1], co);
            if (swap_own[5]) bck_inc_own.n[0] = std::max(bck_inc_own.n[0], co);
            else             bck_inc_own.n[1] = std::max(bck_inc_own.n[1], co);
            stk.push(LFT);
            stk.push(FRO);
            stk.push(BCK);
        }

        Hexa_adap_compute(own, adap_incr);
        assert(adap_incr.cell_adap_incr[own].n[0] == ao);
        assert(adap_incr.cell_adap_incr[own].n[1] == bo);
        assert(adap_incr.cell_adap_incr[own].n[2] == co);
        assert(!is_HX27(own, adap_incr));

        // promote nei faces
        BOT = pNei[0]-1;
        TOP = pNei[1]-1;
        LFT = pNei[2]-1;
        RGT = pNei[3]-1;
        FRO = pNei[4]-1;
        BCK = pNei[5]-1;
        auto& bot_inc_nei = adap_incr.face_adap_incr[BOT];
        auto& top_inc_nei = adap_incr.face_adap_incr[TOP];
        auto& lft_inc_nei = adap_incr.face_adap_incr[LFT];
        auto& rgt_inc_nei = adap_incr.face_adap_incr[RGT];
        auto& fro_inc_nei = adap_incr.face_adap_incr[FRO];
        auto& bck_inc_nei = adap_incr.face_adap_incr[BCK];

        const auto& swap_nei = _canon_info[nei];

        an = adap_incr.cell_adap_incr[nei].n[0];
        bn = adap_incr.cell_adap_incr[nei].n[1];
        cn = adap_incr.cell_adap_incr[nei].n[2];

        if (nei_mod_a) {
            if (swap_nei[0]) bot_inc_nei.n[1] = std::max(bot_inc_nei.n[1], an);
            else             bot_inc_nei.n[0] = std::max(bot_inc_nei.n[0], an);
            if (swap_nei[1]) top_inc_nei.n[1] = std::max(top_inc_nei.n[1], an);
            else             top_inc_nei.n[0] = std::max(top_inc_nei.n[0], an);
            if (swap_nei[4]) fro_inc_nei.n[1] = std::max(fro_inc_nei.n[1], an);
            else             fro_inc_nei.n[0] = std::max(fro_inc_nei.n[0], an);
            if (swap_nei[5]) bck_inc_nei.n[1] = std::max(bck_inc_nei.n[1], an);
            else             bck_inc_nei.n[0] = std::max(bck_inc_nei.n[0], an);
            stk.push(BOT);
            stk.push(TOP);
            stk.push(FRO);
            stk.push(BCK);
        }

        if (nei_mod_b) {
            if (swap_nei[0]) bot_inc_nei.n[0] = std::max(bot_inc_nei.n[0], bn);
            else             bot_inc_nei.n[1] = std::max(bot_inc_nei.n[1], bn);
            if (swap_nei[1]) top_inc_nei.n[0] = std::max(top_inc_nei.n[0], bn);
            else             top_inc_nei.n[1] = std::max(top_inc_nei.n[1], bn);
            if (swap_nei[2]) lft_inc_nei.n[1] = std::max(lft_inc_nei.n[1], bn);
            else             lft_inc_nei.n[0] = std::max(lft_inc_nei.n[0], bn);
            if (swap_nei[3]) rgt_inc_nei.n[1] = std::max(rgt_inc_nei.n[1], bn);
            else             rgt_inc_nei.n[0] = std::max(rgt_inc_nei.n[0], bn);
            stk.push(BOT);
            stk.push(TOP);
            stk.push(RGT);
        }

        if (nei_mod_c) {
            if (swap_nei[2]) lft_inc_nei.n[0] = std::max(lft_inc_nei.n[0], cn);
            else             lft_inc_nei.n[1] = std::max(lft_inc_nei.n[1], cn);
            if (swap_nei[3]) rgt_inc_nei.n[0] = std::max(rgt_inc_nei.n[0], cn);
            else             rgt_inc_nei.n[1] = std::max(rgt_inc_nei.n[1], cn);
            if (swap_nei[4]) fro_inc_nei.n[0] = std::max(fro_inc_nei.n[0], cn);
            else             fro_inc_nei.n[1] = std::max(fro_inc_nei.n[1], cn);
            if (swap_nei[5]) bck_inc_nei.n[0] = std::max(bck_inc_nei.n[0], cn);
            else             bck_inc_nei.n[1] = std::max(bck_inc_nei.n[1], cn);
            stk.push(RGT);
            stk.push(FRO);
            stk.push(BCK);
        }

        Hexa_adap_compute(nei, adap_incr);
        assert(adap_incr.cell_adap_incr[nei].n[0] == an);
        assert(adap_incr.cell_adap_incr[nei].n[1] == bn);
        assert(adap_incr.cell_adap_incr[nei].n[2] == cn);
        assert(!is_HX27(nei, adap_incr));


        // resolve XY problems
        if (!is_FXY(PGi, adap_incr)) continue;

        ao = adap_incr.cell_adap_incr[own].n[0];
        bo = adap_incr.cell_adap_incr[own].n[1];
        co = adap_incr.cell_adap_incr[own].n[2];

        if (is_HX18(own, adap_incr)) {
            if (ao && bo) {

            } else if (ao && co) {

            }
        } else if (is_HX12(own, adap_incr)) {
            if (bo) {

            } else if (co) {

            }
        }

        an = adap_incr.cell_adap_incr[nei].n[0];
        bn = adap_incr.cell_adap_incr[nei].n[1];
        cn = adap_incr.cell_adap_incr[nei].n[2];

        if (is_HX18(nei, adap_incr)) {
            if (an && bn) {

            } else if (an && cn) {

            }
        } else if (is_HX12(nei, adap_incr)) {
            if (bn) {

            } else if (cn) {
                
            }
        }
        

        

    } // end while

    // reorder HX12 in gamma to alpha
    for (E_Int PHi = 0; PHi < hmesh._ng.PHs.size(); PHi++) {
        if (is_HX27(PHi, adap_incr)) {
            std::cout << "HX27?\n";
        }
        if (!hmesh._PHtree.is_enabled(PHi)) continue;
        if (!is_HX12(PHi, adap_incr)) continue;
        E_Int gamma = adap_incr.cell_adap_incr[PHi].n[2];
        if (gamma == 0) {
            K_MESH::Hexahedron::reorder_pgs(hmesh._ng, hmesh._F2E, PHi);
            compute_canon_info_bottom(PHi, 0);
            continue;
        }

        E_Int *pF = hmesh._ng.PHs.get_facets_ptr(PHi);

        K_CONNECT::IdTool::right_shift<6>(pF, 2);
        K_MESH::Hexahedron::reorder_pgs(hmesh._ng, hmesh._F2E, PHi);

        compute_canon_info_bottom(PHi, 0);
        Hexa_adap_compute(PHi, adap_incr);
        assert(is_HX12(PHi, adap_incr) && adap_incr.cell_adap_incr[PHi].n[2] == 0);
    }

    // reorder HX18
    for (E_Int PHi = 0; PHi < hmesh._ng.PHs.size(); PHi++) {
        if (!hmesh._PHtree.is_enabled(PHi)) continue;
        if (!is_HX18(PHi, adap_incr)) continue;

        E_Int *pF = hmesh._ng.PHs.get_facets_ptr(PHi);

        E_Int i0 = 0;
        while (!is_FXY(pF[0]-1, adap_incr)) {
            K_CONNECT::IdTool::right_shift<6>(pF, 1);
            i0++;
            if (i0 == 6) {
                std::cout << "oups\n";
                break;
            }
        }
            
        K_MESH::Hexahedron::reorder_pgs(hmesh._ng, hmesh._F2E, PHi);

        compute_canon_info_bottom(PHi, 0);
        Hexa_adap_compute(PHi, adap_incr);
        assert(is_HX18(PHi, adap_incr) && adap_incr.cell_adap_incr[PHi].n[2] == 0);
    }
} // end function
*/


template <typename mesh_t>
void metric_sensor<mesh_t>::metric_fix(mesh_t& hmesh, output_t& adap_incr)
{
    using cell_incr_t = typename output_t::cell_incr_t;
    using face_incr_t = typename output_t::face_incr_t;

    cell_incr_t c_incr_1, c_incr_0;
    E_Int count;
    E_Int *pF = 0;
    E_Int alpha, beta, gamma;
    std::stack<E_Int> stk;
    E_Int PGi, PHi;

    for (PHi = 0; PHi < hmesh._ng.PHs.size(); PHi++) {
        if (adap_incr.cell_adap_incr[PHi] == 0) continue;
        stk.push(PHi);
    }

    while (!stk.empty()) {
        PHi = stk.top();
        stk.pop();

        count = 0;
        for (E_Int i = 0; i < 3; i++) {
            if (adap_incr.cell_adap_incr[PHi].n[i])
                count++;
        }

        if (count == 0 || count == 3) continue;

        alpha = adap_incr.cell_adap_incr[PHi].n[0];
        beta = adap_incr.cell_adap_incr[PHi].n[1];
        gamma = adap_incr.cell_adap_incr[PHi].n[2];

        const auto& swap = _canon_info[PHi];
        pF = hmesh._ng.PHs.get_facets_ptr(PHi);

        for (int i = 0; i < 6; i++) {
            PGi = pF[i]-1;
            if (!is_FXY(PGi, adap_incr)) continue;

            E_Int nei = NEIGHBOR(PHi, hmesh._F2E, PGi);
            bool ok_nei = (nei != IDX_NONE) && hmesh._PHtree.is_enabled(PHi);

            auto& PG_incr = adap_incr.face_adap_incr[PGi];

            if (count == 1) {
                if (alpha) {
                    if (i != 2 && i != 3) {
                        if (swap[i]) PG_incr.n[0] = 0;
                        else PG_incr.n[1] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    }
                } else if (beta) {
                    if (i == 0 || i == 1) {
                        if (swap[i]) PG_incr.n[1] = 0;
                        else PG_incr.n[0] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    } else if (i == 2 || i == 3) {
                        if (swap[i]) PG_incr.n[0] = 0;
                        else PG_incr.n[1] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    }
                } else if (gamma) {
                    if (i != 0 && i != 1) {
                        if (swap[i]) PG_incr.n[1] = 0;
                        else PG_incr.n[0] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    }
                }
            } else if (count == 2) {
                if (alpha && beta) {
                    if (i != 0 && i != 1) {
                        if (swap[i]) PG_incr.n[0] = 0;
                        else PG_incr.n[1] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    }
                } else if (beta && gamma) {
                    if (i != 2 && i != 3) {
                        if (swap[i]) PG_incr.n[1] = 0;
                        else PG_incr.n[0] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    }
                } else if (gamma && alpha) {
                    if (i == 0 || i == 1) {
                        if (swap[i]) PG_incr.n[0] = 0;
                        else PG_incr.n[1] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    } else if (i == 2 || i == 3) {
                        if (swap[i]) PG_incr.n[1] = 0;
                        else PG_incr.n[0] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    }
                }
            }
        }
    }

    // reorder HX12 in gamma to alpha
    for (PHi = 0; PHi < hmesh._ng.PHs.size(); PHi++) {
        if (!hmesh._PHtree.is_enabled(PHi)) continue;
        if (!is_HX12(PHi, adap_incr)) continue;
        gamma = adap_incr.cell_adap_incr[PHi].n[2];
        if (gamma == 0) continue;

        pF = hmesh._ng.PHs.get_facets_ptr(PHi);

        K_CONNECT::IdTool::right_shift<6>(pF, 2);
        K_MESH::Hexahedron::reorder_pgs(hmesh._ng, hmesh._F2E, PHi);

        compute_canon_info_bottom(PHi, 0);
        Hexa_adap_compute(PHi, adap_incr);
        assert(is_HX12(PHi, adap_incr) && adap_incr.cell_adap_incr[PHi].n[2] == 0);
    }

    // reorder HX18
    for (PHi = 0; PHi < hmesh._ng.PHs.size(); PHi++) {
        if (!hmesh._PHtree.is_enabled(PHi)) continue;
        if (!is_HX18(PHi, adap_incr)) continue;

        pF = hmesh._ng.PHs.get_facets_ptr(PHi);

        E_Int i0 = 0;
        while (!is_FXY(pF[0]-1, adap_incr)) {
            K_CONNECT::IdTool::right_shift<6>(pF, 1);
            i0++;
            if (i0 == 6) {
                std::cout << "oups\n";
                break;
            }
        }
            
        K_MESH::Hexahedron::reorder_pgs(hmesh._ng, hmesh._F2E, PHi);

        compute_canon_info_bottom(PHi, 0);
        Hexa_adap_compute(PHi, adap_incr);
        assert(is_HX18(PHi, adap_incr) && adap_incr.cell_adap_incr[PHi].n[2] == 0);
    }

#ifdef DEBUG_METRIC_SENSOR
    std::cout << "adap_incr fixed\n";
#endif
}

///
template <typename mesh_t>
void metric_sensor<mesh_t>::Q4_adap_compute(E_Int PGi, output_t& adap_incr)
{
  //int h[4];
  E_Float h[4];
  E_Int n1, n2;
  auto& metric = parent_t::_data;
  const E_Int *pN = parent_t::_hmesh._ng.PGs.get_facets_ptr(PGi);

  for (int i = 0; i < 4; i++) {
    n1 = pN[i]; n2 = pN[(i+1)%4];
    h[i] = metric.lengthEval(n1-1, metric[n1-1], n2-1, metric[n2-1]);
  }

  auto f1 = log2(std::min(h[0], h[2]));
  auto f2 = log2(std::min(h[1], h[3]));

  adap_incr.face_adap_incr[PGi].n[0] = round(f1);
  adap_incr.face_adap_incr[PGi].n[1] = round(f2);
}

template <typename mesh_t>
void metric_sensor<mesh_t>::Hexa_adap_compute(E_Int PHi, output_t& adap_incr)
{
  if (PHi == IDX_NONE) return;
  if (!parent_t::_hmesh._PHtree.is_enabled(PHi)) return;

  using face_incr_t = typename output_t::face_incr_t;

  const auto& ng = parent_t::_hmesh._ng;
  const auto& F2E = parent_t::_hmesh._F2E;

  face_incr_t f_incr[6];
  E_Int reorient, i0;
  E_Int PGi;
  const E_Int *pF = ng.PHs.get_facets_ptr(PHi);
  const auto& swap = _canon_info[PHi];

  for (int i = 0; i < 6; i++) {
    PGi = pF[i] - 1;
    if (swap[i] > 0) {
      f_incr[i].n[0] = adap_incr.face_adap_incr[PGi].n[1];
      f_incr[i].n[1] = adap_incr.face_adap_incr[PGi].n[0];
    } else {
      f_incr[i].n[0] = adap_incr.face_adap_incr[PGi].n[0];
      f_incr[i].n[1] = adap_incr.face_adap_incr[PGi].n[1];
    }
  }

  using cell_incr_t = typename output_t::cell_incr_t;

  cell_incr_t cell_incr(0);

  E_Int incrX[4] = {f_incr[0].n[0], f_incr[4].n[0], f_incr[1].n[0], f_incr[5].n[0]};
  cell_incr.n[0] = std::max(0,*std::min_element(incrX, incrX + 4));

  E_Int incrY[4] = {f_incr[0].n[1], f_incr[2].n[0], f_incr[1].n[1], f_incr[3].n[0]};
  cell_incr.n[1] = std::max(0,*std::min_element(incrY, incrY + 4));

  E_Int incrZ[4] = {f_incr[2].n[1], f_incr[4].n[1], f_incr[3].n[1], f_incr[5].n[1]};
  cell_incr.n[2] = std::max(0,*std::min_element(incrZ, incrZ + 4));

  adap_incr.cell_adap_incr[PHi] = cell_incr;
}

///
template <typename mesh_t>
bool metric_sensor<mesh_t>::fill_adap_incr(output_t& adap_incr, bool do_agglo)
{
  using cell_incr_t = typename output_t::cell_incr_t; //int_tuple<3>
  using face_incr_t = typename output_t::face_incr_t; // int_tuple<2>
  
  // init face_adap_incr
  adap_incr.face_adap_incr.clear();
  E_Int nb_faces = parent_t::_hmesh._ng.PGs.size();
  adap_incr.face_adap_incr.resize(nb_faces, face_incr_t(0));

  // init cell_adap_incr
  adap_incr.cell_adap_incr.clear();
  E_Int nb_cells = parent_t::_hmesh._ng.PHs.size();
  adap_incr.cell_adap_incr.resize(nb_cells, cell_incr_t(0));

  _single_pass_done = true; // to exit after first outer loop (based on sensor update)

  E_Int PHi, PGi, stride;
  
  std::vector<E_Int> hx27;

  for (PHi = 0; PHi < nb_cells; PHi++) {
    if (!parent_t::_hmesh._PHtree.is_enabled(PHi)) continue;
    stride = parent_t::_hmesh._ng.PHs.stride(PHi);
    const E_Int *pF = parent_t::_hmesh._ng.PHs.get_facets_ptr(PHi);
    for (int j = 0; j < stride; j++) {
      PGi = pF[j] - 1;
      // only hexa for now
      Q4_adap_compute(PGi, adap_incr);
    }
    compute_canon_info_bottom(PHi, 0);
    Hexa_adap_compute(PHi, adap_incr);
  }

#ifdef DEBUG_METRIC_SENSOR
  std::cout << "adap_incr filled\n";
#endif

  metric_fix(parent_t::_hmesh, adap_incr);

  E_Int cmax = 0;
  for (size_t k = 0; (k < adap_incr.cell_adap_incr.size()); ++k)
    cmax = std::max(cmax, adap_incr.cmax(k));

  E_Int cmin = IDX_NONE; // max int32
  for (size_t k = 0; (k < adap_incr.cell_adap_incr.size()) && (cmin != 0); ++k)
    cmin = std::min(cmin, adap_incr.cmin(k));

  //std::cout << "cell adapt incr min/max : " << cmin << "/" << cmax << std::endl;

  return (cmin != 0 || cmax != 0);
}

}
#endif
