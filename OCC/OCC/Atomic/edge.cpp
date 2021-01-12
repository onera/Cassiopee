// ########## Builder and destructors for planar objects ###############
static void edge_dealloc(py_edge* self)
{
    self->ob_type->tp_free((PyObject*)self);
}
// Create a new edge
static PyObject* edge_new(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
    py_edge* self = (py_edge*)type->tp_alloc(type, 0);
    return (PyObject*)self;
}

// ############################################################################
// Initializer for the Plane incident wave instances :
// ---------------------------------------------------
static const char edge_init_doc[] =
    u8R"RAW(
        Edge creation.
)RAW";
static int edge_init(py_edge* self, PyObject* args)
{
    return 0;
}
// ############################################################################
// Accessors and modifiers
// -----------------------
static char edge_get_n_doc[] = R"RAW(
...
)RAW";
static PyObject* edge_get_n(py_edge* self, void* closure)
{
    //return Py_BuildValue( "d", double( self->fwave->n_wave ) );
    Py_RETURN_NONE;
}
static int edge_set_n_wave(py_edge* self, PyObject* value, void* closure)
{
    return 0;
}
// ############################################################################
// Getter and setter table (attribs):
// ----------------------------------
static PyGetSetDef edge_gettersetters[] = {
    //{(char*)"n", (getter)edge_get_n, (setter)edge_set_n, edge_doc, NULL},
    {NULL}};

// ############################################################################
// method table :
// -------------
static PyMethodDef edge_methods[] = {
    //{(char*)"pressure", (PyCFunction)py_planar_pressure, METH_VARARGS, planar_pressure_doc},
    {NULL}};
// ############################################################################
static PyObject* edge_repr(py_edge* self)
{
    return PyString_FromString("Object representation");
}
// print method
static PyObject* edge_str(py_edge* self)
{
    return PyString_FromString("Object content");
}
// ############################################################################
static const char planar_doc[] = R"RAW(
Class for CAD edge.
)RAW";
#if defined(_WIN64) || defined(_WIN32)
__declspec( dllexport )
#endif
PyTypeObject EdgeType = {
        PyObject_HEAD_INIT( NULL ) 0,                 // ob_size
        "Occ.edge",                                   // nom Python du type
        sizeof(py_edge),                              // Taille basique de l'objet
        0,                                            // Taille de l'"item"
        (destructor)edge_dealloc,                     // Destructeur
        (printfunc)0,                                 // Pour afficher l'objet
        0,                                            // getattr
        0,                                            // setattr
        0,                                            // compare
        (reprfunc)edge_repr,                          // Representation de l'objet
        0,                                            // L'objet considere comme un nombre
        0,                                            // L'objet considere comme une séquence
        0,                                            // L'objet considere comme un dictionnaire
        0,                                            // L'objet considere comme une table de hashage
        0,                                            // L'objet considere comme une fonction
        (reprfunc)edge_str,                           // Représentation chaîne caractere de l'objet
        0,                                            // get_attro
        0,                                            // set_attro
        0,                                            // L'objet comme buffer
        Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,     // Objet banal ici ;)
        edge_doc,                                     // Documentation de la classe
        0,                                            // tp_traverse
        0,                                            // tp_clear
        0,                                            /* tp_richcompare */
        0,                                            /* tp_weaklistoffset */
        0,                                            /* tp_iter */
        0,                                            /* tp_iternext */
        edge_methods,                                 /* tp_methods */
        0,                                            /* tp_members */
        edge_gettersetters,                           /* tp_getset */
        0,                                            /* tp_base */
        0,                                            /* tp_dict */
        0,                                            /* tp_descr_get */
        0,                                            /* tp_descr_set */
        0,                                            /* tp_dictoffset */
        (initproc)edge_init,                          /* tp_init */
        0,                                            /* tp_alloc */
        edge_new,                                     /* tp_new */
};
