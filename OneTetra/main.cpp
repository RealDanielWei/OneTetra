#include "Tetra_system.h"


int main() {
	Tetra_system::Tetra_domain tdomain("nodes_cube.txt", "elements_cube.txt");
	cout << "output Sh...Ne=" << tdomain.get_Ne() << endl;
	tdomain.output_Sh();
	cout << "output Se...Nh=" << tdomain.get_Nh() << endl;
	tdomain.output_Se();
}