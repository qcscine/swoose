/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_AMINOACIDDATAFORTESTS_H
#define SWOOSE_AMINOACIDDATAFORTESTS_H

#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>

namespace Scine {
namespace Swoose {

namespace AminoAcidDataForTests {

std::string ALA = "HETATM    1  N   ALA     0      -0.966   0.493   1.500  0.00  0.00           N \n"
                  "HETATM    2  CA  ALA     0       0.257   0.418   0.692  0.00  0.00           C \n"
                  "HETATM    3  C   ALA     0      -0.094   0.017  -0.716  0.00  0.00           C \n"
                  "HETATM    4  O   ALA     0      -1.056  -0.682  -0.923  0.00  0.00           O \n"
                  "HETATM    5  CB  ALA     0       1.204  -0.620   1.296  0.00  0.00           C \n"
                  "HETATM    6  OXT ALA     0       0.661   0.439  -1.742  0.00  0.00           O \n"
                  "END ";
std::string ARG = "HETATM    1  N   ARG     0      -0.469   1.110  -0.993  0.00  0.00           N \n"
                  "HETATM    2  CA  ARG     0       0.004   2.294  -1.708  0.00  0.00           C \n"
                  "HETATM    3  C   ARG     0      -0.907   2.521  -2.901  0.00  0.00           C \n"
                  "HETATM    4  O   ARG     0      -1.827   1.789  -3.242  0.00  0.00           O \n"
                  "HETATM    5  CB  ARG     0       1.475   2.150  -2.127  0.00  0.00           C \n"
                  "HETATM    6  CG  ARG     0       1.745   1.017  -3.130  0.00  0.00           C \n"
                  "HETATM    7  CD  ARG     0       3.210   0.954  -3.557  0.00  0.00           C \n"
                  "HETATM    8  NE  ARG     0       4.071   0.726  -2.421  0.00  0.00           N \n"
                  "HETATM    9  CZ  ARG     0       5.469   0.624  -2.528  0.00  0.00           C \n"
                  "HETATM   10  NH1 ARG     0       6.259   0.404  -1.405  0.00  0.00           N \n"
                  "HETATM   11  NH2 ARG     0       6.078   0.744  -3.773  0.00  0.00           N1+ \n"
                  "HETATM   12  OXT ARG     0      -0.588   3.659  -3.574  0.00  0.00           O \n"
                  "END ";
std::string ASN = "HETATM    1  N   ASN     0      -0.293   1.686   0.094  0.00  0.00           N\n"
                  "HETATM    2  CA  ASN     0      -0.448   0.292  -0.340  0.00  0.00           C  \n"
                  "HETATM    3  C   ASN     0      -1.846  -0.179  -0.031  0.00  0.00           C  \n"
                  "HETATM    4  O   ASN     0      -2.510   0.402   0.794  0.00  0.00           O  \n"
                  "HETATM    5  CB  ASN     0       0.562  -0.588   0.401  0.00  0.00           C  \n"
                  "HETATM    6  CG  ASN     0       1.960  -0.197  -0.002  0.00  0.00           C  \n"
                  "HETATM    7  OD1 ASN     0       2.132   0.697  -0.804  0.00  0.00           O  \n"
                  "HETATM    8  ND2 ASN     0       3.019  -0.841   0.527  0.00  0.00           N  \n"
                  "HETATM    9  OXT ASN     0      -2.353  -1.243  -0.673  0.00  0.00           O  \n"
                  "END";
std::string ASP = "HETATM    1  N   ASP     0      -0.317   1.688   0.066  0.00  0.00           N \n"
                  "HETATM    2  CA  ASP     0      -0.470   0.286  -0.344  0.00  0.00           C \n"
                  "HETATM    3  C   ASP     0      -1.868  -0.180  -0.029  0.00  0.00           C \n"
                  "HETATM    4  O   ASP     0      -2.534   0.415   0.786  0.00  0.00           O \n"
                  "HETATM    5  CB  ASP     0       0.539  -0.580   0.413  0.00  0.00           C \n"
                  "HETATM    6  CG  ASP     0       1.938  -0.195   0.004  0.00  0.00           C \n"
                  "HETATM    7  OD1 ASP     0       2.109   0.681  -0.810  0.00  0.00           O \n"
                  "HETATM    8  OD2 ASP     0       2.992  -0.826   0.543  0.00  0.00           O \n"
                  "HETATM    9  OXT ASP     0      -2.374  -1.256  -0.652  0.00  0.00           O \n"
                  "END ";

std::string CYS = "HETATM    1  N   CYS     0       1.585   0.483  -0.081  0.00  0.00           N \n"
                  "HETATM    2  CA  CYS     0       0.141   0.450   0.186  0.00  0.00           C \n"
                  "HETATM    3  C   CYS     0      -0.095   0.006   1.606  0.00  0.00           C \n"
                  "HETATM    4  O   CYS     0       0.685  -0.742   2.143  0.00  0.00           O \n"
                  "HETATM    5  CB  CYS     0      -0.533  -0.530  -0.774  0.00  0.00           C \n"
                  "HETATM    6  SG  CYS     0      -0.247   0.004  -2.484  0.00  0.00           S \n"
                  "HETATM    7  OXT CYS     0      -1.174   0.443   2.275  0.00  0.00           O \n"
                  "END ";

std::string GLN = "HETATM    1  N   GLN     0       1.858  -0.148   1.125  0.00  0.00           N \n"
                  "HETATM    2  CA  GLN     0       0.517   0.451   1.112  0.00  0.00           C \n"
                  "HETATM    3  C   GLN     0      -0.236   0.022   2.344  0.00  0.00           C \n"
                  "HETATM    4  O   GLN     0      -0.005  -1.049   2.851  0.00  0.00           O \n"
                  "HETATM    5  CB  GLN     0      -0.236  -0.013  -0.135  0.00  0.00           C \n"
                  "HETATM    6  CG  GLN     0       0.529   0.421  -1.385  0.00  0.00           C \n"
                  "HETATM    7  CD  GLN     0      -0.213  -0.036  -2.614  0.00  0.00           C \n"
                  "HETATM    8  OE1 GLN     0      -1.252  -0.650  -2.500  0.00  0.00           O \n"
                  "HETATM    9  NE2 GLN     0       0.277   0.236  -3.839  0.00  0.00           N \n"
                  "HETATM   10  OXT GLN     0      -1.165   0.831   2.878  0.00  0.00           O \n"
                  "END ";

std::string GLU = "HETATM    1  N   GLU     0       1.199   1.867  -0.117  0.00  0.00           N \n"
                  "HETATM    2  CA  GLU     0       1.138   0.515   0.453  0.00  0.00           C \n"
                  "HETATM    3  C   GLU     0       2.364  -0.260   0.041  0.00  0.00           C \n"
                  "HETATM    4  O   GLU     0       3.010   0.096  -0.916  0.00  0.00           O \n"
                  "HETATM    5  CB  GLU     0      -0.113  -0.200  -0.062  0.00  0.00           C \n"
                  "HETATM    6  CG  GLU     0      -1.360   0.517   0.461  0.00  0.00           C \n"
                  "HETATM    7  CD  GLU     0      -2.593  -0.187  -0.046  0.00  0.00           C \n"
                  "HETATM    8  OE1 GLU     0      -2.485  -1.161  -0.753  0.00  0.00           O \n"
                  "HETATM    9  OE2 GLU     0      -3.811   0.269   0.287  0.00  0.00           O \n"
                  "HETATM   10  OXT GLU     0       2.737  -1.345   0.737  0.00  0.00           O \n"
                  "END ";

std::string GLY = "HETATM    1  N   GLY     0       1.931   0.090  -0.034  0.00  0.00           N \n"
                  "HETATM    2  CA  GLY     0       0.761  -0.799  -0.008  0.00  0.00           C \n"
                  "HETATM    3  C   GLY     0      -0.498   0.029  -0.005  0.00  0.00           C \n"
                  "HETATM    4  O   GLY     0      -0.429   1.235  -0.023  0.00  0.00           O \n"
                  "HETATM    5  OXT GLY     0      -1.697  -0.574   0.018  0.00  0.00           O \n"
                  "END ";

std::string HIS = "HETATM    1  N   HIS     0      -0.040  -1.210   0.053  0.00  0.00           N \n"
                  "HETATM    2  CA  HIS     0       1.172  -1.709   0.652  0.00  0.00           C \n"
                  "HETATM    3  C   HIS     0       1.083  -3.207   0.905  0.00  0.00           C \n"
                  "HETATM    4  O   HIS     0       0.040  -3.770   1.222  0.00  0.00           O \n"
                  "HETATM    5  CB  HIS     0       1.484  -0.975   1.962  0.00  0.00           C \n"
                  "HETATM    6  CG  HIS     0       2.940  -1.060   2.353  0.00  0.00           C \n"
                  "HETATM    7  ND1 HIS     0       3.380  -2.075   3.129  0.00  0.00           N1+ \n"
                  "HETATM    8  CD2 HIS     0       3.960  -0.251   2.046  0.00  0.00           C \n"
                  "HETATM    9  CE1 HIS     0       4.693  -1.908   3.317  0.00  0.00           C \n"
                  "HETATM   10  NE2 HIS     0       5.058  -0.801   2.662  0.00  0.00           N \n"
                  "HETATM   11  OXT HIS     0       2.247  -3.882   0.744  0.00  0.00           O \n"
                  "END ";

std::string ILE = "HETATM    1  N   ILE     0      -1.944   0.335  -0.343  0.00  0.00           N \n"
                  "HETATM    2  CA  ILE     0      -0.487   0.519  -0.369  0.00  0.00           C \n"
                  "HETATM    3  C   ILE     0       0.066  -0.032  -1.657  0.00  0.00           C \n"
                  "HETATM    4  O   ILE     0      -0.484  -0.958  -2.203  0.00  0.00           O \n"
                  "HETATM    5  CB  ILE     0       0.140  -0.219   0.814  0.00  0.00           C \n"
                  "HETATM    6  CG1 ILE     0      -0.421   0.341   2.122  0.00  0.00           C \n"
                  "HETATM    7  CG2 ILE     0       1.658  -0.027   0.788  0.00  0.00           C \n"
                  "HETATM    8  CD1 ILE     0       0.206  -0.397   3.305  0.00  0.00           C \n"
                  "HETATM    9  OXT ILE     0       1.171   0.504  -2.197  0.00  0.00           O \n"
                  "END ";

std::string LEU = "HETATM    1  N   LEU     0      -1.661   0.627  -0.406  0.00  0.00           N \n"
                  "HETATM    2  CA  LEU     0      -0.205   0.441  -0.467  0.00  0.00           C \n"
                  "HETATM    3  C   LEU     0       0.180  -0.055  -1.836  0.00  0.00           C \n"
                  "HETATM    4  O   LEU     0      -0.591  -0.731  -2.474  0.00  0.00           O \n"
                  "HETATM    5  CB  LEU     0       0.221  -0.583   0.585  0.00  0.00           C \n"
                  "HETATM    6  CG  LEU     0      -0.170  -0.079   1.976  0.00  0.00           C \n"
                  "HETATM    7  CD1 LEU     0       0.256  -1.104   3.029  0.00  0.00           C \n"
                  "HETATM    8  CD2 LEU     0       0.526   1.254   2.250  0.00  0.00           C \n"
                  "HETATM    9  OXT LEU     0       1.382   0.254  -2.348  0.00  0.00           O \n"
                  "END ";

std::string LYS = "HETATM    1  N   LYS     0       1.422   1.796   0.198  0.00  0.00           N \n"
                  "HETATM    2  CA  LYS     0       1.394   0.355   0.484  0.00  0.00           C \n"
                  "HETATM    3  C   LYS     0       2.657  -0.284  -0.032  0.00  0.00           C \n"
                  "HETATM    4  O   LYS     0       3.316   0.275  -0.876  0.00  0.00           O \n"
                  "HETATM    5  CB  LYS     0       0.184  -0.278  -0.206  0.00  0.00           C \n"
                  "HETATM    6  CG  LYS     0      -1.102   0.282   0.407  0.00  0.00           C \n"
                  "HETATM    7  CD  LYS     0      -2.313  -0.351  -0.283  0.00  0.00           C \n"
                  "HETATM    8  CE  LYS     0      -3.598   0.208   0.329  0.00  0.00           C \n"
                  "HETATM    9  NZ  LYS     0      -4.761  -0.400  -0.332  0.00  0.00           N1+ \n"
                  "HETATM   10  OXT LYS     0       3.050  -1.476   0.446  0.00  0.00           O \n"
                  "END ";

std::string MET = "HETATM    1  N   MET     0      -1.816   0.142  -1.166  0.00  0.00           N \n"
                  "HETATM    2  CA  MET     0      -0.392   0.499  -1.214  0.00  0.00           C \n"
                  "HETATM    3  C   MET     0       0.206   0.002  -2.504  0.00  0.00           C \n"
                  "HETATM    4  O   MET     0      -0.236  -0.989  -3.033  0.00  0.00           O \n"
                  "HETATM    5  CB  MET     0       0.334  -0.145  -0.032  0.00  0.00           C \n"
                  "HETATM    6  CG  MET     0      -0.273   0.359   1.277  0.00  0.00           C \n"
                  "HETATM    7  SD  MET     0       0.589  -0.405   2.678  0.00  0.00           S \n"
                  "HETATM    8  CE  MET     0      -0.314   0.353   4.056  0.00  0.00           C \n"
                  "HETATM    9  OXT MET     0       1.232   0.661  -3.066  0.00  0.00           O \n"
                  "END ";

std::string PHE = "HETATM    1  N   PHE     0       1.317   0.962   1.014  0.00  0.00           N \n"
                  "HETATM    2  CA  PHE     0      -0.020   0.426   1.300  0.00  0.00           C \n"
                  "HETATM    3  C   PHE     0      -0.109   0.047   2.756  0.00  0.00           C \n"
                  "HETATM    4  O   PHE     0       0.879  -0.317   3.346  0.00  0.00           O \n"
                  "HETATM    5  CB  PHE     0      -0.270  -0.809   0.434  0.00  0.00           C \n"
                  "HETATM    6  CG  PHE     0      -0.181  -0.430  -1.020  0.00  0.00           C \n"
                  "HETATM    7  CD1 PHE     0       1.031  -0.498  -1.680  0.00  0.00           C \n"
                  "HETATM    8  CD2 PHE     0      -1.314  -0.018  -1.698  0.00  0.00           C \n"
                  "HETATM    9  CE1 PHE     0       1.112  -0.150  -3.015  0.00  0.00           C \n"
                  "HETATM   10  CE2 PHE     0      -1.231   0.333  -3.032  0.00  0.00           C \n"
                  "HETATM   11  CZ  PHE     0      -0.018   0.265  -3.691  0.00  0.00           C \n"
                  "HETATM   12  OXT PHE     0      -1.286   0.113   3.396  0.00  0.00           O \n"
                  "END ";

std::string PRO = "HETATM    1  N   PRO     0      -0.816   1.108   0.254  0.00  0.00           N \n"
                  "HETATM    2  CA  PRO     0       0.001  -0.107   0.509  0.00  0.00           C \n"
                  "HETATM    3  C   PRO     0       1.408   0.091   0.005  0.00  0.00           C \n"
                  "HETATM    4  O   PRO     0       1.650   0.980  -0.777  0.00  0.00           O \n"
                  "HETATM    5  CB  PRO     0      -0.703  -1.227  -0.286  0.00  0.00           C \n"
                  "HETATM    6  CG  PRO     0      -2.163  -0.753  -0.439  0.00  0.00           C \n"
                  "HETATM    7  CD  PRO     0      -2.218   0.614   0.276  0.00  0.00           C \n"
                  "HETATM    8  OXT PRO     0       2.391  -0.721   0.424  0.00  0.00           O \n"
                  "END ";

std::string SER = "HETATM    1  N   SER     0       1.525   0.493  -0.608  0.00  0.00           N \n"
                  "HETATM    2  CA  SER     0       0.100   0.469  -0.252  0.00  0.00           C \n"
                  "HETATM    3  C   SER     0      -0.053   0.004   1.173  0.00  0.00           C \n"
                  "HETATM    4  O   SER     0       0.751  -0.760   1.649  0.00  0.00           O \n"
                  "HETATM    5  CB  SER     0      -0.642  -0.489  -1.184  0.00  0.00           C \n"
                  "HETATM    6  OG  SER     0      -0.496  -0.049  -2.535  0.00  0.00           O \n"
                  "HETATM    7  OXT SER     0      -1.084   0.440   1.913  0.00  0.00           O \n"
                  "END ";

std::string THR = "HETATM    1  N   THR     0       1.543  -0.702   0.430  0.00  0.00           N \n"
                  "HETATM    2  CA  THR     0       0.122  -0.706   0.056  0.00  0.00           C \n"
                  "HETATM    3  C   THR     0      -0.038  -0.090  -1.309  0.00  0.00           C \n"
                  "HETATM    4  O   THR     0       0.732   0.761  -1.683  0.00  0.00           O \n"
                  "HETATM    5  CB  THR     0      -0.675   0.104   1.079  0.00  0.00           C \n"
                  "HETATM    6  OG1 THR     0      -0.193   1.448   1.103  0.00  0.00           O \n"
                  "HETATM    7  CG2 THR     0      -0.511  -0.521   2.466  0.00  0.00           C \n"
                  "HETATM    8  OXT THR     0      -1.039  -0.488  -2.110  0.00  0.00           O \n"
                  "END ";

std::string TRP = "HETATM    1  N   TRP     0       1.278   1.121   2.059  0.00  0.00           N \n"
                  "HETATM    2  CA  TRP     0      -0.008   0.417   1.970  0.00  0.00           C \n"
                  "HETATM    3  C   TRP     0      -0.490   0.076   3.357  0.00  0.00           C \n"
                  "HETATM    4  O   TRP     0       0.308  -0.130   4.240  0.00  0.00           O \n"
                  "HETATM    5  CB  TRP     0       0.168  -0.868   1.161  0.00  0.00           C \n"
                  "HETATM    6  CG  TRP     0       0.650  -0.526  -0.225  0.00  0.00           C \n"
                  "HETATM    7  CD1 TRP     0       1.928  -0.418  -0.622  0.00  0.00           C \n"
                  "HETATM    8  CD2 TRP     0      -0.186  -0.256  -1.396  0.00  0.00           C \n"
                  "HETATM    9  NE1 TRP     0       1.978  -0.095  -1.951  0.00  0.00           N \n"
                  "HETATM   10  CE2 TRP     0       0.701   0.014  -2.454  0.00  0.00           C \n"
                  "HETATM   11  CE3 TRP     0      -1.564  -0.210  -1.615  0.00  0.00           C \n"
                  "HETATM   12  CZ2 TRP     0       0.190   0.314  -3.712  0.00  0.00           C \n"
                  "HETATM   13  CZ3 TRP     0      -2.044   0.086  -2.859  0.00  0.00           C \n"
                  "HETATM   14  CH2 TRP     0      -1.173   0.348  -3.907  0.00  0.00           C \n"
                  "HETATM   15  OXT TRP     0      -1.806   0.001   3.610  0.00  0.00           O \n"
                  "END ";

std::string TYR = "HETATM    1  N   TYR     0       1.320   0.952   1.428  0.00  0.00           N \n"
                  "HETATM    2  CA  TYR     0      -0.018   0.429   1.734  0.00  0.00           C \n"
                  "HETATM    3  C   TYR     0      -0.103   0.094   3.201  0.00  0.00           C \n"
                  "HETATM    4  O   TYR     0       0.886  -0.254   3.799  0.00  0.00           O \n"
                  "HETATM    5  CB  TYR     0      -0.274  -0.831   0.907  0.00  0.00           C \n"
                  "HETATM    6  CG  TYR     0      -0.189  -0.496  -0.559  0.00  0.00           C \n"
                  "HETATM    7  CD1 TYR     0       1.022  -0.589  -1.219  0.00  0.00           C \n"
                  "HETATM    8  CD2 TYR     0      -1.324  -0.102  -1.244  0.00  0.00           C \n"
                  "HETATM    9  CE1 TYR     0       1.103  -0.282  -2.563  0.00  0.00           C \n"
                  "HETATM   10  CE2 TYR     0      -1.247   0.210  -2.587  0.00  0.00           C \n"
                  "HETATM   11  CZ  TYR     0      -0.032   0.118  -3.252  0.00  0.00           C \n"
                  "HETATM   12  OH  TYR     0       0.044   0.420  -4.574  0.00  0.00           O \n"
                  "HETATM   13  OXT TYR     0      -1.279   0.184   3.842  0.00  0.00           O \n"
                  "END ";

std::string VAL = "HETATM    1  N   VAL     0       1.564  -0.642   0.454  0.00  0.00           N \n"
                  "HETATM    2  CA  VAL     0       0.145  -0.698   0.079  0.00  0.00           C \n"
                  "HETATM    3  C   VAL     0      -0.037  -0.093  -1.288  0.00  0.00           C \n"
                  "HETATM    4  O   VAL     0       0.703   0.784  -1.664  0.00  0.00           O \n"
                  "HETATM    5  CB  VAL     0      -0.682   0.086   1.098  0.00  0.00           C \n"
                  "HETATM    6  CG1 VAL     0      -0.497  -0.528   2.487  0.00  0.00           C \n"
                  "HETATM    7  CG2 VAL     0      -0.218   1.543   1.119  0.00  0.00           C \n"
                  "HETATM    8  OXT VAL     0      -1.022  -0.529  -2.089  0.00  0.00           O \n"
                  "END ";

inline std::pair<std::string, int> getAminoAcid(const std::string& residueName) {
  if (residueName == "ALA")
    return std::make_pair(ALA, 13);
  else if (residueName == "ARG")
    return std::make_pair(ARG, 27);
  else if (residueName == "ASN")
    return std::make_pair(ASN, 17);
  else if (residueName == "ASP")
    return std::make_pair(ASP, 16);
  else if (residueName == "CYS")
    return std::make_pair(CYS, 14);
  else if (residueName == "GLN")
    return std::make_pair(GLN, 20);
  else if (residueName == "GLU")
    return std::make_pair(GLU, 19);
  else if (residueName == "GLY")
    return std::make_pair(GLY, 10);
  else if (residueName == "HIS")
    return std::make_pair(HIS, 20);
  else if (residueName == "ILE")
    return std::make_pair(ILE, 22);
  else if (residueName == "LEU")
    return std::make_pair(LEU, 22);
  else if (residueName == "LYS")
    return std::make_pair(LYS, 24);
  else if (residueName == "MET")
    return std::make_pair(MET, 20);
  else if (residueName == "PHE")
    return std::make_pair(PHE, 23);
  else if (residueName == "PRO")
    return std::make_pair(PRO, 17);
  else if (residueName == "SER")
    return std::make_pair(SER, 14);
  else if (residueName == "THR")
    return std::make_pair(THR, 17);
  else if (residueName == "TRP")
    return std::make_pair(TRP, 27);
  else if (residueName == "TYR")
    return std::make_pair(TYR, 24);
  else if (residueName == "VAL")
    return std::make_pair(VAL, 19);
  else
    throw std::runtime_error("Unsupported Amino Acid");
}

} // namespace AminoAcidDataForTests
} // namespace Swoose
} // namespace Scine

#endif // SWOOSE_REFERENCEDATAFORTESTS_H
