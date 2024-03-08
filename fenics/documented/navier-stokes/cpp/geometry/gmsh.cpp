#include <gmsh.h>
#include <set>
#include <iostream>


int main(int argc, char **argv){
    gmsh::initialize();

    gmsh::model::add("aboba");
    double lc = 1e-1;

    int pt1 = gmsh::model::occ::addPoint(0, 0, 0, lc);
    int pt2 = gmsh::model::occ::addPoint(1, 0, 0, lc);
    int pt3 = gmsh::model::occ::addPoint(1, 1, 0, lc);
    int pt4 = gmsh::model::occ::addPoint(0, 1, 0, lc);

    int ln1 = gmsh::model::occ::addLine(pt1, pt2);
    int ln2 = gmsh::model::occ::addLine(pt2, pt3);
    int ln3 = gmsh::model::occ::addLine(pt3, pt4);
    int ln4 = gmsh::model::occ::addLine(pt4, pt1);

    int surf = gmsh::model::occ::addCurveLoop({ln1, ln2, ln3, ln4});
    int circ = gmsh::model::occ::addCircle(0.5, 0.5, 0, 0.25, -1);
    int circ_loop = gmsh::model::occ::addCurveLoop({circ});
    int plane = gmsh::model::occ::addPlaneSurface({surf, circ_loop});

    int phys = gmsh::model::addPhysicalGroup(2, {plane});

    gmsh::model::occ::synchronize();

    gmsh::model::mesh::generate(2);
    
   // gmsh::option::setNumber("Mesh.SaveAll", 1);
    gmsh::write("../mesh/kek.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();
    gmsh::finalize();
    return 0;
}