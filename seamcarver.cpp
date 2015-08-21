#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <utility>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <set>

#include "CImg.h"

double luminosity(int r, int g, int b);
double pythagorean(double a, double b);
bool valid_coord(std::pair<int, int> coord, int width, int height);
std::vector<std::vector<double> > calculate_luminosity(cimg_library::CImg<unsigned char> &image);
std::vector<std::vector<std::pair<double, double> > > calculate_gradient(std::vector<std::vector<double> > luminosity_map);
double calculate_seam(std::vector<std::vector<std::pair<double, double> > > &gradient_map, int x, int y);
std::vector<std::pair<int, int> > trace_seam(std::vector<std::vector<std::pair<double, double> > > gradient_map, int index);

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cout << "Not enough arguments, Usage: ./seamcarver <path_to_image> <number_of_seams_to_carve>" << "\n";
        return 1;
    }

    std::cout << "Loading image" << "\n";
    std::string fname = argv[1];
    cimg_library::CImg<unsigned char> image(fname.c_str());
    cimg_library::CImg<unsigned char> draw_image(fname.c_str()); //TODO figure out if i can undo the draw without having second copy

    std::cout << "Calculating luminosity" << "\n";
    std::vector<std::vector<double> > luminosity_map = calculate_luminosity(image);

    std::vector<std::vector<int> > energy_map(luminosity_map.size(), std::vector<int>(luminosity_map[0].size(), 0));
    cimg_library::CImgDisplay energy(draw_image, "Manually set special energy");

    char red[3] = { 255, 0, 0 };
    char green[3] = { 0, 255, 0 };
    int cursor_size = 5;
    while (!energy.is_closed()) {
        energy.wait();
        if (energy.button()) {
            std::vector<char> current_color;
            int energy_val = 0;
            if (energy.button() & 1) { // left click
                current_color = std::vector<char>(red, red + 3);
                energy_val = -1;
            } else if (energy.button() & 2) { // right click
                current_color = std::vector<char>(green, green + 3);
                energy_val = 1;
            } else if (energy.button() & 3) { // middle click
                //todo figure out if you can uncolor image
                //current_color = ?
                energy_val = 0; 
            }
            for (int i = -cursor_size; i <= cursor_size; i++) {
                for (int j = -cursor_size; j <= cursor_size; j++) {
                    energy_map[energy.mouse_y() + i][energy.mouse_x() + j] = energy_val;
                    draw_image.draw_point(energy.mouse_x() + j, energy.mouse_y() + i, &current_color[0]);
                }
            }
        }
        draw_image.display(energy);
    } 

    std::vector<std::vector<std::vector<char> > > im_vec;
    for (int i = 0; i < image.height(); i++) {
        std::vector<std::vector<char> > row;
        for (int j = 0; j < image.width(); j++) {
            std::vector<char> color;
            color.push_back(image.atXY(j, i, 0, 0));
            color.push_back(image.atXY(j, i, 0, 1));
            color.push_back(image.atXY(j, i, 0, 2));
            row.push_back(color);
        }
        im_vec.push_back(row);
    }

    int numseam = atoi(argv[2]);
    for (int master = 0; master < numseam; master++) {
        std::cout << "Computing and removing seam #" << master + 1 << "\r" << std::flush;

        std::vector<std::vector<std::pair<double, double> > > gradient_map = calculate_gradient(luminosity_map);

        //set all energies specified in energy_map to negative/positive on the gradient map
        for (int i = 0; i < energy_map.size(); i++) {
            for (int j = 0; j < energy_map[0].size(); j++) {
                if (energy_map[i][j] != 0) {
                    gradient_map[i][j].second = energy_map[i][j] * 1000;
                }
            }
        }

        double min_seam = DBL_MAX;
        int loc = -1;
        for (int i = 0; i < gradient_map.back().size(); i++) {
            calculate_seam(gradient_map, gradient_map.size() - 1, i);
            if (gradient_map.back()[i].first < min_seam) {
                min_seam = gradient_map.back()[i].first;
                loc = i;
            }
        }

        //remove seams and the object manipulation map coordinates
        std::vector<std::pair<int, int> > toremove = trace_seam(gradient_map, loc);
        for (int i = 0; i < toremove.size(); i++) {
            int row = toremove[i].first;
            int col = toremove[i].second;
            im_vec[row].erase(im_vec[row].begin() + col);
            luminosity_map[row].erase(luminosity_map[row].begin() + col);
            energy_map[row].erase(energy_map[row].begin() + col);
        }
    }

    cimg_library::CImg<unsigned char> newimage(im_vec[0].size(), im_vec.size(), 1, 3);
    for (int i = 0; i < im_vec.size(); i++) {
        for (int j = 0; j < im_vec[i].size(); j++) {
            char color[3] = { im_vec[i][j][0], im_vec[i][j][1], im_vec[i][j][2] };
            newimage.draw_point(j, i, color);
        }
    }

    std::cout << "\n" << "Computation finished, displaying image" << "\n"; 

    cimg_library::CImgDisplay display(newimage, "Seam Carved Image");

    while (!display.is_closed()) {
        display.wait();

        newimage.display(display);
    }
}

//pair<double, double> --> current energy seam, gradient value
double calculate_seam(std::vector<std::vector<std::pair<double, double> > > &gradient_map, int x, int y) {
    if (x == 0) return gradient_map[x][y].second;
    if (gradient_map[x][y].first < DBL_MAX) return gradient_map[x][y].first;

    std::vector<double> res(3, DBL_MAX);

    for (int i = -1; i <= 1; i++) {
        if (valid_coord(std::pair<int, int>(x - 1, y + i), gradient_map[0].size(), gradient_map.size())) {
            res[i + 1] = calculate_seam(gradient_map, x - 1, y + i);
        }
    }

    gradient_map[x][y].first = gradient_map[x][y].second + std::min(std::min(res[0], res[1]), res[2]);

    return gradient_map[x][y].first;
}

//returns list of coordinates of pixels to remove
std::vector<std::pair<int, int> > trace_seam(std::vector<std::vector<std::pair<double, double> > > gradient_map, int index) {
    std::vector<std::pair<int, int> > toremove;
    for (int i = gradient_map.size() - 1; i >= 0; i--) {
        toremove.push_back(std::pair<int, int>(i, index));
        int loc = -1;
        double val = DBL_MAX;
        for (int j = -1; j <= 1; j++) {
            if (valid_coord(std::pair<double, double>(i - 1, index + j), gradient_map[0].size(), gradient_map.size())) {
                if (gradient_map[i - 1][index + j].first < val) {
                    val = gradient_map[i - 1][index + j].first;
                    loc = index + j;
                }
            }
        } 
        index = loc;
    }
    return toremove;
}

double luminosity(int r, int g, int b) {
    return  0.2126 * r + 0.7152 * g + 0.0722 * b;
}

std::vector<std::vector<double> > calculate_luminosity(cimg_library::CImg<unsigned char> &image) {
    std::vector<std::vector<double> > luminosity_map(image.height(), std::vector<double>(image.width(), DBL_MAX));
    for (int i = 0; i < image.height(); i++) {
        for (int j = 0; j < image.width(); j++) {
            luminosity_map[i][j] = luminosity(image.atXY(j, i, 0, 0), image.atXY(j, i, 0, 1), image.atXY(j, i, 0, 2));
        }
    }
    return luminosity_map;
}

//pair<double, double> --> current energy seam, gradient value
std::vector<std::vector<std::pair<double, double> > > calculate_gradient(std::vector<std::vector<double> > luminosity_map) {
    std::vector<std::vector<std::pair<double, double> > > gradient_map(luminosity_map.size(), 
                                                                       std::vector<std::pair<double, double> >(luminosity_map[0].size(),
                                                                                           std::pair<double, double>(DBL_MAX, DBL_MAX)));
    int width = luminosity_map[0].size();
    int height = luminosity_map.size();
    for (int i = 0; i < luminosity_map.size(); i++) {
        for (int j = 0; j < luminosity_map[i].size(); j++) {
            double current_val = luminosity_map[i][j];

            double u = valid_coord(std::pair<int, int>(i - 1, j), width, height) ? luminosity_map[i - 1][j] : current_val;
            double d = valid_coord(std::pair<int, int>(i + 1, j), width, height) ? luminosity_map[i + 1][j] : current_val;
            double l = valid_coord(std::pair<int, int>(i, j - 1), width, height) ? luminosity_map[i][j - 1] : current_val;
            double r = valid_coord(std::pair<int, int>(i, j + 1), width, height) ? luminosity_map[i][j + 1] : current_val;

            double dx = r - l;
            double dy = u - d;
            gradient_map[i][j] = std::pair<double, double>(DBL_MAX, pythagorean(dx, dy));
            if (i == 0) {
                gradient_map[i][j].first = gradient_map[i][j].second;
            }
        }
    }
    return gradient_map;
}

double pythagorean(double a, double b) {
    return sqrt(a * a + b * b);
}

bool valid_coord(std::pair<int, int> coord, int width, int height) {
    if (coord.first < 0 || coord.second < 0) return false;
    return coord.first < height && coord.second < width;
}
