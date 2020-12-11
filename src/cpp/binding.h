#ifndef BINDING_H_INCLUDED
#define BINDING_H_INCLUDED


extern "C" {
    OrthOpt* OrthOpt() {return new OrthOpt(); }
    OrthOpt* OrthOpt(Mesh* m) {return new OrthOpt(m); }

    Mesh* Mesh() {return new Mesh(); }

    void load_vertices_array(Mesh mesh, double vertices[], size_t s) {
        Mesh::load_elements_array(vertices, s);
    }

    void load_elements_array(Mesh mesh, double vertices[], size_t s) {
        Mesh::load_elements_array(vertices, s);
    }

    unsigned int optimize(OrthOpt* obj, double penalization_power,
                          unsigned int max_it) {

        return code
    }

}

#endif // BINDING_H_INCLUDED
