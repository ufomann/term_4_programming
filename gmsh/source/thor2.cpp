#include <set>
#include <gmsh.h>
#include <math.h>
#include <iostream>

#define dim first
#define tag second

#define PI 3.1415926

int draw_circ(double R, int n_on_arc, double r, double lc){
    double ang = 2 * PI / n_on_arc;

    int center = gmsh::model::occ::addPoint(R, 0, 0, lc);

    std::vector<int> pts(n_on_arc);
    for (int i = 0; i < n_on_arc; i++){
        pts[i] = gmsh::model::occ::addPoint((r * cos(i * ang) + R), 0,
                                    r * sin(i * ang), lc);
    }

    std::vector<int> arcs(n_on_arc);
    for (int i = 0; i < n_on_arc; i++){
        arcs[i] = gmsh::model::occ::addCircleArc(pts[i], center, pts[(i + 1) % n_on_arc]);
    }

    gmsh::model::occ::remove({{0, center}});

    int circ = gmsh::model::occ::addCurveLoop(arcs);
    return circ;
}

int main(int argc, char **argv){
    double r_in = 1.0;
    double r_out = 1.5;
    double R = 5.0;
    double lc = 0.1;
    int n_on_arc = 10;
    double ang = 2 * PI / n_on_arc;
    double thor_ang = 2 * PI;
    double n_rings = 100;

    gmsh::initialize();

    gmsh::model::add("thor2");

    int center = gmsh::model::occ::addPoint(R, 0, 0, lc);

    int in_circ = draw_circ(R, 20, r_in, lc);
    int out_circ = draw_circ(R, 20, r_out, lc);
    
    int slice = gmsh::model::occ::addPlaneSurface({in_circ, out_circ});
    
    gmsh::vectorpair out;
    gmsh::model::occ::revolve({{2, slice}}, 0, 0, 0, 0, 0, 1, thor_ang, out, {n_rings});

    gmsh::vectorpair s;
    gmsh::model::getEntities(s, 2);
    std::vector<int> sl;
    for(auto surf : s) sl.push_back(surf.second);
    int l = gmsh::model::geo::addSurfaceLoop(sl);
    gmsh::model::geo::addVolume({l});


    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(3);
 
    gmsh::write("mesh/thor2.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();
    gmsh::finalize();
    return 0;
}