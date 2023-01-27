#include "Mesh.h"
#include "IO.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <iomanip> //set precision

//TODO: look at vertice id order when declaring element


void IO::load_mesh_auto(std::string f_mesh, 
                        std::string f_vertices, 
                        std::string f_elements) {
    if (f_vertices.size() != 0) {
        std::cout << "Read tetgen format mesh" << std::endl;
        load_vertices_tetgen(f_vertices);
        load_elements_tetgen(f_elements);
        return;
    }
    
    std::string ext = f_mesh.substr(f_mesh.find_last_of(".")+1);
    if ( (ext.compare("mesh") == 0) or (ext.compare("meshb") == 0) ) {
        std::cout << "Read Medit mesh format" << std::endl;
        IO::load_mesh_medit(f_mesh);
    }
    else if (ext.compare("ugi") == 0) {
        std::cout << "Read PFLOTRAN ugi mesh format" << std::endl;
        IO::load_mesh_PFLOTRAN(f_mesh);
    }
    else if (ext.compare("dat") == 0) {
        std::cout << "Read Salome DAT mesh format" << std::endl;
        IO::load_mesh_DAT_salome(f_mesh);
    }
    else {
        std::cerr << std::endl;
        std::cerr << "Read Mesh - ERROR" << std::endl;
        std::cerr << "Extension \"" << ext << "\" not recognized" << std::endl;
    }
    return;
}

void IO::save_mesh_auto(std::string f_out) {
    std::cout << "Output file: " << f_out << std::endl;
    std::string ext = f_out.substr(f_out.find_last_of(".")+1);
    if (ext.compare("xyz") == 0) {
        std::cout << "Save optimized vertices in xyz format" << std::endl;
        IO::save_vertices_xyz(f_out);
    }
    else if ( (ext.compare("mesh") == 0) or (ext.compare("meshb") == 0) ) {
        std::cout << "Save optimized mesh in Medit format" << std::endl;
        IO::save_mesh_medit(f_out);
    }
    else if (ext.compare("dat") == 0) {
        std::cout << "Save optimized mesh in Salome DAT format" << std::endl;
        IO::save_mesh_DAT_salome(f_out);
    }
    else if (ext.compare("ugi") == 0) {
        std::cout << "Save optimized mesh in PFLOTRAN ugi format" << std::endl;
        IO::save_mesh_PFLOTRAN(f_out);
    }
    else {
        std::cerr << std::endl;
        std::cerr << "Write Mesh - ERROR" << std::endl;
        std::cerr << "Extension \"" << ext << "\" not recognized" << std::endl;
        std::cerr << "Write mesh in xyz format instead" << std::endl;
        IO::save_vertices_tetgen(f_out);
    }
    return;
}


// XYZ //
void IO::load_vertices_xyz(std::string filename) {
    // X Y Z file
    double x, y, z;
    unsigned int vertice_id = 0;
    std::ifstream src(filename);
    while (src >> x >> y >> z) {
        vertice_id += 1;
        mesh->add_vertice(x,y,z,vertice_id);
    }
    src.close();
}
void IO::save_vertices_xyz(std::string filename) {
    std::ofstream out(filename);
    out << std::scientific;
    out << std::setprecision(8);
    for (Vertice* v : mesh->vertices) {
        //for (size_t i=0; i<mesh->dim; i++) 
        out << (Eigen::RowVectorXd) *(v->coor) << std::endl;
    }
    out.close();
}
void IO::load_elements_with_vertice_ids(std::string filename) {
    //each line one element, and list of id separated by a space delimiter
    std::ifstream src(filename);
    std::string delimiter = " ";
    size_t pos = 0;
    std::string token;
    for( std::string line; std::getline( src, line ); ) {
        //one line one element
        std::vector<unsigned int> ids;
        while ((pos = line.find(delimiter)) != std::string::npos) {
            token = line.substr(0, pos);
            line.erase(0, pos + delimiter.length());
            ids.push_back(std::stoi(token));
        }
        mesh->add_element(ids);
    }
    src.close();
}


// TetGen //
void IO::load_vertices_tetgen(std::string filename) {
    // X Y Z file
    double x, y, z;
    unsigned int vertice_id;
    std::ifstream src(filename);
    while (src >> vertice_id >> x >> y >> z) {
        mesh->add_vertice(x,y,z,vertice_id);
    }
    src.close();
}
void IO::save_vertices_tetgen(std::string filename) {
    std::ofstream src(filename);
    unsigned int count = 1;
    for (Vertice* v : mesh->vertices) {
        src << std::fixed;
        src << count << ' ';
        src << std::scientific;
        src << std::setprecision(8);
        src << (v->coor)[0] << ' ';
        src << (v->coor)[1] << ' ';
        src << (v->coor)[2] << std::endl;
        count++;
    }
    src.close();
}
void IO::load_elements_tetgen(std::string filename) {
    //each line one element, and list of id separated by a space delimiter
    std::ifstream src(filename);
    unsigned int read_value;
    std::vector<unsigned int> ids;
    while (src >> read_value) {
        ids.push_back(read_value);
        if (ids.size() == 5) {
            ids.erase(ids.begin());
            mesh->add_element(ids);
            ids.clear();
        }
    }
    src.close();
}


// Medit //
void IO::load_mesh_medit(std::string filename) {
    double x, y, z;
    unsigned int vertice_id = 1;
    std::ifstream src(filename);
    std::string keyword = "NONE";
    unsigned int line_num = 0;
    std::vector<unsigned int> element_ids;
    std::string token;
    std::istringstream f;
    unsigned int val;

    for ( std::string text; std::getline( src, text ); ) {
        line_num++;
        if ( s_len_trim ( text ) == 0 ) {
            keyword = "NONE";
            continue;
        }
        if ( text[0] == '#' ) {
            continue;
        }
        
        if ( s_eqi ( text, "CORNERS" ) ) {
            keyword = "SKIP";
        }
        else if ( s_begin ( text, "DIMENSION 2" ) ) {
            mesh->dim = 2;
        }
        else if ( s_begin ( text, "DIMENSION 3" ) ) {
            mesh->dim = 3;
        }
        else if ( s_eqi ( text, "EDGES" ) ) {
            keyword = "SKIP";
        }
        else if ( s_eqi ( text, "END" ) ) {
            break;
        }
        else if ( s_begin ( text, "MESHVERSIONFORMATTED" ) ) {}
        else if ( s_eqi ( text, "NORMALATQUADRILATERALVERTICES" ) ) {
            keyword = "SKIP";
        }
        else if ( s_eqi ( text, "NORMALATTRIANGLEVERTICES" ) ) {
            keyword = "SKIP";
        }
        else if ( s_eqi ( text, "NORMALATVERTICES" ) ) {
            keyword = "SKIP";
        }
        else if ( s_eqi ( text, "NORMALS" ) ) {
            keyword = "SKIP";
        }
        else if ( s_eqi ( text, "REQUIREDEDGES" ) ) {
            keyword = "SKIP";
        }
        else if ( s_eqi ( text, "REQUIREDVERTICES" ) ) {
            keyword = "SKIP";
        }
        else if ( s_eqi ( text, "RIDGES" ) ) {
            keyword = "SKIP";
        }
        else if ( s_eqi ( text, "TANGENTATEDGES" ) ) {
            keyword = "SKIP";
        }
        else if ( s_eqi ( text, "TANGENTS" ) ) {
            keyword = "SKIP";
        }
        else if ( s_eqi ( text, "TETRAHEDRA" ) ||
                  s_eqi ( text, "TETRAHEDRONS" ) ) {
            if (mesh->dim == 3) keyword = "TETRAHEDRONS";
            else keyword = "SKIP";
        }
        else if ( s_eqi ( text, "PRISMS" ) ) {
            if (mesh->dim == 3) keyword = "PRISMS";
            else keyword = "SKIP";
        }
        else if ( s_eqi ( text, "PYRAMIDS" ) ) {
            if (mesh->dim == 3) keyword = "PYRAMIDS";
            else keyword = "SKIP";
        }
        else if ( s_eqi ( text, "HEXAHEDRA" ) ||
                  s_eqi ( text, "HEXAHEDRONS" ) ) {
            if (mesh->dim == 3) keyword = "HEXAHEDRONS";
            else keyword = "SKIP";
        }
        else if ( s_eqi ( text, "TRIANGLES" ) ) {
            if (mesh->dim == 2) keyword = "TRIANGLES";
            else keyword = "SKIP";
        }
        else if ( s_eqi ( text, "QUADRILATERALS" ) ) {
            if (mesh->dim == 2) keyword = "QUADRILATERALS";
            else keyword = "SKIP";
        }
        else if ( s_eqi ( text, "VERTICES" ) ) {
            keyword = "VERTICES";
        }
        //
        //  Presumably, numeric data to be processed by keyword.
        //
        else if ( s_eqi ( keyword, "QUADRILATERALS" ) )
        {
            atoi ( text.c_str() ); //number of quad
            keyword = "QUADRILATERAL_VERTEX";
        }
        else if ( s_eqi ( keyword, "QUADRILATERAL_VERTEX" ) )
        {
            element_ids.clear();
            f.str(text);
            for (int n=0; n<4; n++) {
                f >> val;
                element_ids.push_back(val);
            }
            mesh->add_element(element_ids);
        }
        else if ( s_eqi ( keyword, "TRIANGLES" ) )
        {
            atoi ( text.c_str() ); //number of triangle
            keyword = "TRIANGLE_VERTEX";
        }
        else if ( s_eqi ( keyword, "TRIANGLE_VERTEX" ) )
        {
            element_ids.clear();
            f.str(text);
            for (int n=0; n<3; n++) {
                f >> val;
                element_ids.push_back(val);
            }
            mesh->add_element(element_ids);
        }
        else if ( s_eqi ( keyword, "HEXAHEDRONS" ) ) {
            atoi ( text.c_str ( ) ); //number of hex
            keyword = "HEXAHEDRON_VERTEX";
        }
        else if ( s_eqi ( keyword, "HEXAHEDRON_VERTEX" ) ) {
            element_ids.clear();
            f.str(text);
            for (int n=0; n<8; n++) {
                f >> val;
                element_ids.push_back(val);
            }
            mesh->add_element(element_ids);
        }
        else if ( s_eqi ( keyword, "TETRAHEDRONS" ) ) {
            atoi ( text.c_str ( ) ); //number of tet
            keyword = "TETRAHEDRON_VERTEX";
        }
        else if ( s_eqi ( keyword, "TETRAHEDRON_VERTEX" ) ) {
            element_ids.clear();
            f.str(text);
            for (int n=0; n<4; n++) {
                f >> val;
                element_ids.push_back(val);
            }
            mesh->add_element(element_ids);
        }
        else if ( s_eqi ( keyword, "PRISMS" ) ) {
            atoi ( text.c_str ( ) ); //number of prisms
            keyword = "PRISMS_VERTEX";
        }
        else if ( s_eqi ( keyword, "PRISMS_VERTEX" ) ) {
            element_ids.clear();
            f.str(text);
            for (int n=0; n<6; n++) {
                f >> val;
                element_ids.push_back(val);
            }
            mesh->add_element(element_ids);
        }
        else if ( s_eqi ( keyword, "PYRAMIDS" ) ) {
            atoi ( text.c_str ( ) ); //number of pyramids
            keyword = "PYRAMIDS_VERTEX";
        }
        else if ( s_eqi ( keyword, "PYRAMIDS_VERTEX" ) ) {
            element_ids.clear();
            f.str(text);
            for (int n=0; n<5; n++) {
                f >> val;
                element_ids.push_back(val);
            }
            mesh->add_element(element_ids);
        }
        else if ( s_eqi ( keyword, "VERTICES" ) ) {
            atoi ( text.c_str ( ) ); //number of vertices
            keyword = "VERTEX_COORDINATE";
        }
        else if ( s_eqi ( keyword, "VERTEX_COORDINATE" ) ) {
            f.str(text);
            if (mesh->dim == 3) {
                f >> x >> y >> z;
                mesh->add_vertice(x,y,z,vertice_id);
            } 
            else {
                f >> x >> y;
                mesh->add_vertice(x,y,vertice_id);
            }
            vertice_id++;
        }
        else if ( s_eqi ( keyword, "SKIP" ) ) {}
        else {
            std::cerr << std::endl;
            std::cerr << "Medit IO - Fatal error!" << std::endl;
            std::cerr << "  Could not find keyword while reading line "
                      << line_num << std::endl;
            std::cerr << "\"" << text << "\"." << std::endl;
            exit ( 1 );
        }
    }
    //
    //  Close the file.
    //
   src.close();

    return;
}

int IO::s_len_trim ( std::string s ) {
    //from https://people.math.sc.edu/Burkardt/cpp_src/medit_io/medit_io.cpp
    int n;
    n = s.length ( );
    while ( 0 < n ) {
        if ( s[n-1] != ' ' && s[n-1] != '\n' ) {
          return n;
        }
        n = n - 1;
    }
    return n;
}

bool IO::s_eqi ( std::string s1, std::string s2 ) {
    //from https://people.math.sc.edu/Burkardt/cpp_src/medit_io/medit_io.cpp
    int i;
    int nchar;
    int s1_length;
    int s2_length;

    s1_length = s1.length ( );
    s2_length = s2.length ( );

    if ( s1_length < s2_length )
    {
      nchar = s1_length;
    }
    else
    {
      nchar = s2_length;
    }
    //
    //  The strings are not equal if they differ over their common length.
    //
    for ( i = 0; i < nchar; i++ ) {
        if ( (char) std::toupper(s1[i]) != (char) std::toupper(s2[i]) ) {
            return false;
        }
    }
    //
    //  The strings are not equal if the longer one includes nonblanks
    //  in the tail.
    //
    if ( nchar < s1_length ) {
        for ( i = nchar; i < s1_length; i++ ) {
            if ( s1[i] != ' ' ) {
                return false;
            }
        }
    }
    else if ( nchar < s2_length ) {
        for ( i = nchar; i < s2_length; i++ ) {
            if ( s2[i] != ' ' ) {
              return false;
            }
        }
    }
    return true;
}

bool IO::s_begin ( std::string s1, std::string s2 ) {
    //from https://people.math.sc.edu/Burkardt/cpp_src/medit_io/medit_io.cpp
    int i;
    int n1;
    int n2;
    n1 = s1.length ( );
    n2 = s2.length ( );
    if ( n1 < n2 ) {
        return false;
    }
    for ( i = 0; i < n2; i++ ) {
        if ( (char) std::toupper(s1[i]) != (char) std::toupper(s2[i]) ) {
            return false;
        }
    }
    return true;
}

void IO::save_mesh_medit(std::string filename) {
    std::ofstream out;
    std::string ext = filename.substr(filename.find_last_of(".")+1);
    if (ext.compare("mesh") == 0) {
        out.open(filename);
    }
    else if (ext.compare("meshb") == 0) {
        out.open(filename, std::ofstream::binary);
    }
    if ( !out )
    {
        std::cerr << std::endl;
        std::cerr << "Medit write - Fatal error!" << std::endl;
        std::cerr << "  Unable to open output file" << std::endl;
        exit ( 1 );
    }
    out << "MeshVersionFormatted 2" << std::endl;
    out << std::endl;
    out << "Dimension " << mesh->dim << std::endl;
    out << std::endl;
    out << "# Created by OrthOpt: " << std::endl;
    out << "# The Mesh Orthogonality Optimizer" << std::endl;
    out << std::endl;
    //out elements
    unsigned int count;
    int elem_types[6] = {-4, -3, 4, 5, 6, 8};
    for (auto elem_type : elem_types) {
        //count elements
        count = 0;
        for (Element* elem : mesh->elements) {
            if (elem->type == elem_type) count++;
        }
        if (count != 0) {
            if (elem_type == -3) {
                out << "Triangles" << std::endl;
            }
            else if (elem_type == -4) {
                out << "Quadrilaterals" << std::endl;
            }
            else if (elem_type == 4) {
                out << "Tetrahedra" << std::endl;
            }
            else if (elem_type == 5) {
                out << "Pyramids" << std::endl;
            }
            else if (elem_type == 6) {
                out << "Prisms" << std::endl;
            }
            else if (elem_type == 8) {
                out << "Hexahedra" << std::endl;
            }
            out << count << std::endl;
            for (Element* elem : mesh->elements) {
                if (elem->type != elem_type) {
                    continue;
                }
                for (unsigned int id : elem->vertice_ids) {
                    out << id << ' ';
                }
                out << "1" << std::endl;
            }
        }
    }
    //out vertices
    out << std::endl;
    out << "Vertices" << std::endl;
    out << mesh->n_vertices << std::endl;
    out << std::scientific;
    out << std::setprecision(8);
    for (Vertice* v : mesh->vertices) {
        out << (Eigen::RowVectorXd) *(v->coor) << " 1" << std::endl;
    }

    out << std::endl;
    out << "End";
    out.close();
    return;
}


// PFLOTRAN ugi format//
void IO::load_mesh_PFLOTRAN(std::string filename) {
    // X Y Z file
    double x, y, z;
    unsigned int n_vertices, n_elements, n_v;
    std::string type;
    unsigned int id=0;
    std::ifstream src(filename);
    unsigned int count = 0;
    unsigned int i = 0;
    std::string line;
    std::istringstream f;
    
    while ( getline(src, line) ) {
        if (count == 0) { 
            f.str(line);
            f >> n_elements;
            f >> n_vertices;
            count++;
            continue;
        }
        if (count < n_elements+1) {
            type = line[0];
            line.erase(0,2);
            if (type.compare("T") == 0) n_v = 4;
            else if (type.compare("P") == 0) n_v = 5;
            else if (type.compare("W") == 0) n_v = 6;
            else if (type.compare("H") == 0) n_v = 8;
            else {
                std::cerr << "Unknown element type: " << type << std::endl;
                exit(1);
            }
            f.clear();
            f.str(line);
            std::vector<unsigned int> ids;
            for (size_t j=0; j<n_v; j++) {
                f >> id;
                ids.push_back(id);
            }
            mesh->add_element(ids);
        }
        else {
            f.clear();
            f.str(line);
            f >> x >> y >> z;
            mesh->add_vertice(x,y,z,i+1);
            i++;
        }
        count++;
    }
    src.close();
}
void IO::save_mesh_PFLOTRAN(std::string filename) {
    std::ofstream out(filename);
    // header
    out << mesh->n_elements << ' ';
    out << mesh->n_vertices << std::endl;
    //elements
    for (Element* e : mesh->elements) {
        if (e->type == 4) out << 'T';
        if (e->type == 5) out << 'P';
        if (e->type == 6) out << 'W';
        if (e->type == 8) out << 'H';
        for (unsigned int id : e->vertice_ids) {
            out << id << ' ';
        }
        out << std::endl;
    }
    //vertices
    out << std::scientific;
    out << std::setprecision(8);
    for (Vertice* v : mesh->vertices) {
        out << (v->coor)[0] << ' ';
        out << (v->coor)[1] << ' ';
        out << (v->coor)[2] << std::endl;
    }
    out.close();
}

void IO::load_mesh_DAT_salome(std::string filename) {
    //Note: in 2D, only X and Y coordinate are read
    std::ifstream src(filename);
    unsigned int n_v, n_e;
    unsigned int id, code, n_nodes, id_node;
    double x=0., y=0., z=0.;
    std::istringstream f;
    std::string line;
    src >> n_v >> n_e;
    for (unsigned int i=0; i<n_v; i++) {
        if (mesh->dim == 2) src >> id >> x >> y;
        else src >> id >> x >> y >> z;
        mesh->add_vertice(x,y,z,id);
    }
    getline(src, line); //flush the rest of the previous line
    for (unsigned int i=0; i<n_e; i++) {
        getline(src, line);
        f.clear();
        f.str(line);
        f >> id >> code;
        if (mesh->dim == 2) {
            switch (code) {
                case 203:
                    n_nodes = 3;
                    break;
                case 204:
                    n_nodes = 4;
                    break;
                default:
                    n_nodes = 0;
                    break;
            }
        }
        else {
            switch (code) {
                case 304:
                    n_nodes = 4;
                    break;
                case 305:
                    n_nodes = 5;
                    break;
                case 306:
                    n_nodes = 6;
                    break;
                case 308:
                    n_nodes = 8;
                    break;
                default:
                    n_nodes = 0;
                    break;
            }
        }
        if (n_nodes) {
            std::vector<unsigned int> ids;
            for (size_t j=0; j<n_nodes; j++) {
                f >> id_node;
                ids.push_back(id_node);
            }
            mesh->add_element(ids);
        }
    }
    src.close();
}
void IO::save_mesh_DAT_salome(std::string filename) {
    std::ofstream out(filename);
    // header
    out << mesh->n_vertices << ' ';
    out << mesh->n_elements << std::endl;
    //vertices
    out << std::scientific;
    out << std::setprecision(8);
    for (Vertice* v : mesh->vertices) {
        out << v->natural_id << ' ';
        out << (v->coor)[0] << ' ';
        out << (v->coor)[1] << ' ';
        out << (v->coor)[2] << std::endl;
    }
    out.close();
    //elements
    for (Element* e : mesh->elements) {
        out << e->natural_id << ' ';
        if (e->type == -3) out << 203;
        if (e->type == -4) out << 204;
        if (e->type == 4) out << 304;
        if (e->type == 5) out << 305;
        if (e->type == 6) out << 306;
        if (e->type == 8) out << 308;
        for (unsigned int id : e->vertice_ids) {
            out << id << ' ';
        }
        out << std::endl;
    }
    out.close();
}


// VTK //

// GMSH //



