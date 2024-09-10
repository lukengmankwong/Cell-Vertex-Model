#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <random>

class Tissue;

#include "libraries.h"
#include "tissue.h"
#include "vertex.h"
#include "edge.h"
#include "cell.h"

double random(double min, double max, unsigned int seed);

void outputData(const Tissue& Tissue);

void writeCellsFile(Tissue* T, const std::string& filename_cells);
void writeDirectorsFile(Tissue* T, const std::string& filename_directors);
void writeCellDefectsFile(Tissue* T, const std::string& filename_cell_defects);
void writeVertexDefectsFile(Tissue* T, const std::string& filename_vertex_defects);

#endif // FUNCTIONS_H
