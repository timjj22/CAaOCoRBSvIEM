/*
  implementation of compression scheme as outlined in "Collision-Aware and Online Compression of Rigid Body Simulations via Integrated Error Minimization"

  Inputs and output are all DMAT objects as outlined in the libigl library.
*/

#include <fstream>
#include <iostream>

#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>

#include "json.hpp"

#include "Compressor.h"
#include "Decompressor.h"

using json = nlohmann::json;

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    std::cout << "specify config.json file to use" << std::endl;
    return 1;
  }

  // using nlohmann's json library
  /*
    JSON keys:

    'ascii_output'                 : bool for the binary/ascii dmat
    'input_file'                   : input simulation dmat
    'compressed_output_dir'        : directory to output the compressed files
    'compressed_output_file'       : file names to output the compressed representation
    'decompressed_output_file'     : file name for the decompressed representation
    'contact_info'                 : optional argument for the contact information
    'decompress_rate'              : playback rate for the decompressed output
    'compression_params'           : parameters for the error bounds
         'position_threshold'
	 'position_peak_threshold'
	 'rotation_threshold'
	 'rotation_peak_threshold'
    'format'                       : list of strings for the ordering of the dataset (p or q)
  */
  json j;
  std::ifstream i(argv[1], std::ifstream::in);
  i >> j;

  // READ IN ALL OF THE SIMULATION DATA AND CONFIGURATIONS

  // output formatting for the DMAT
  bool ascii = j["ascii_output"];

  // get the input/output files
  std::string input  = j["input_file"];
  std::string output_dir = j["compressed_output_dir"];
  std::string output_c = j["compressed_output_file"];
  std::string output_d = j["decompressed_output_file"];

  // get the format for the numbers in the input file ['q', 'p'] denotes having the quaternion and then position per object
  std::vector<std::string> formatting = j["format"];

  // get the optional contact info
  std::string contact_info = "";
  if(j.count("contact_info") == 1)
    contact_info = j["contact_info"];

  // use the compression parameters to construct the Compressor. If 1 is specified, all should be
  double position_threshold = 0, position_peak_threshold = 0, rotation_threshold = 0, rotation_peak_threshold = 0;
  if(j["compression_params"].count("position_threshold"))
  {
    position_threshold      = j["compression_params"]["position_threshold"];
    position_peak_threshold = j["compression_params"]["position_peak_threshold"];
    rotation_threshold      = j["compression_params"]["rotation_threshold"];
    rotation_peak_threshold = j["compression_params"]["rotation_peak_threshold"];
  }
  else
    std::cout << "Error parsing the threshold parameters" << std::endl;

  // now read in the input file line-by-line and compress
  Eigen::MatrixXd simData;
  if(!igl::readDMAT(input, simData))
  {
    std::cout << "Error reading input file" << std::endl;
    return 1;
  }
  
  // 7 elements per object (4 for quat, 3 for pos), but first row is the timesteps for the simulation
  int32_t objects = (simData.rows() - 1) / 7;

  // read in the optional compression information
  // provide a custom comparator to deal with floating point keys
  auto cmp = [](double a, double b){ return (std::abs(a - b) > 1e-7) && (a < b); };
  std::map<double, std::vector<std::pair<int32_t, int32_t> >, decltype(cmp)> contact_pairs(cmp);
  if(contact_info != "")
  {
    std::ifstream cInfo(contact_info, std::ifstream::in);
    std::string line;
    while(std::getline(cInfo, line))
    {
      // now want to split up for the commas
      double time;
      int32_t o1, o2;
      char c1, c2;
      std::stringstream ss(line);
      ss >> time >> c1 >> o1 >> c2 >> o2; // parse the row
      contact_pairs[time].emplace_back(o1, o2);
    }
  }


  // CREATE COMPRESSION OBJECTS (ONE PER SIMULATION OBJECT)
  

  // make a compressor/decompressor pair per object
  std::vector<Compressor> comps;
  std::vector<Decompressor> decomps;
  for(int32_t i = 0; i < objects; i++)
  {
    comps.emplace_back(position_threshold, position_peak_threshold, rotation_threshold, rotation_peak_threshold);
    decomps.emplace_back();
  }


  // COMPRESS THE DATA FRAME BY FRAME (THIS CAN ALSO BE SETUP FOR STREAMING)


  // this is trivially parallelizable, reading in the column corresponding to a single simulation frame, across the rows of each object
  /*
    Want to record a given contact if an object goes into or out of persistent contact
    state     :  0  0  1  1  1  0  0
    state - 1 :     0  0  1  1  1  0
    tripped   :  0  0  1  0  0  1  0
    
    trip the contact when: state XOR (state - 1) == 1
  */
  Eigen::Matrix<bool, Eigen::Dynamic, 1> contactStatus = Eigen::Matrix<bool, Eigen::Dynamic, 1>::Zero(objects); // 1 represents this object is in contact
  Eigen::Matrix<bool, Eigen::Dynamic, 1> prevContactStatus = Eigen::Matrix<bool, Eigen::Dynamic, 1>::Zero(objects); // 1 represents this object is in contact
  for(int32_t i = 0; i < simData.cols(); i++)
  {
    // for the new time, update the contactStatus
    contactStatus *= 0;
    if(contact_pairs.count(simData(0, i)))
    {
      for(auto k: contact_pairs[simData(0, i)])
      {
	contactStatus(k.first) = 1;
	contactStatus(k.second) = 1;
      }
    }
    for(int32_t j = 0; j < objects; j++)
    {
      if(contactStatus(j) ^ prevContactStatus(j))
	std::cout << "contact_info" << std::endl;
      
      comps[j].compressFrame(CompressionEngine::Frame(simData(0, i),
						      formatting[0] == "q" ? simData.col(i).segment(1 + 7*(j), 4) : simData.col(i).segment(1 + 7*(j) + 3, 4),
						      formatting[0] == "p" ? simData.col(i).segment(1 + 7*(j), 3) : simData.col(i).segment(1 + 7*(j) + 4, 3)
						      ),
			     // 0, // if we have contact information, use this
			     (contactStatus(j) ^ prevContactStatus(j)),
			     (i < simData.cols() - 1) ? 0 : 1 // force the final frame to be added into the compressed representation
			     );
    }
    
    prevContactStatus = contactStatus;
  }


  // READ OUT THE COMPRESSED DATA

  
  // now that the data has been read in, output the compressed representation per object (kinda messy)
  std::vector<Eigen::MatrixXd> compressedPositions;
  std::vector<Eigen::MatrixXd> compressedRotations;
  for(int32_t i = 0; i < objects; i++)
  {
    compressedPositions.push_back(comps[i].getCompressedMatrixPositionRepresentation());
    compressedRotations.push_back(comps[i].getCompressedMatrixRotationRepresentation());
    igl::writeDMAT(output_dir + std::to_string(i) + "P" + output_c, comps[i].getCompressedMatrixPositionRepresentation(), ascii);
    igl::writeDMAT(output_dir + std::to_string(i) + "Q" + output_c, comps[i].getCompressedMatrixRotationRepresentation(), ascii);
  }


  // PASS THE COMPRESSED DATA TO THE DECOMPRESSOR

  
  // populate the decompressor data
  for(int32_t i = 0; i < objects; i++)
  {
    // set the data to point at for the decompressor
    decomps[i].setCompressedData(&compressedRotations[i], &compressedPositions[i]);
  }

  // decompress the inputs, loop through the time bounds at a fixed timestep
  double decompressRate = j["decompress_rate"];
  double endTime = simData(0, simData.cols() - 1);
  double startTime = simData(0, 0);

  std::vector<std::vector<CompressionEngine::Frame>> decompressedFrame; // this is confusing, but the outer level is the timestep, inner level is each object

  int32_t timeSteps = (endTime - startTime) / decompressRate;
  Eigen::MatrixXd states(1 + 7*objects, timeSteps); 


  // DECOMPRESS THE DATA

  
  for(int32_t i = 0; i < timeSteps; i++)
  {
    decompressedFrame.push_back(std::vector<CompressionEngine::Frame>());
    // loop through each object
    for(int32_t j = 0; j < objects; j++)
    {
      CompressionEngine::Frame objectFrame = decomps[j].stateAtTime(i * decompressRate);
      decompressedFrame.back().push_back(objectFrame);

      // put it into a large state matrix, same format as the input
      states(0, i) = objectFrame.time;
      if(formatting[0] == "q")
      {
	states.col(i).segment(1 + 7*(j) + 0, 4) = objectFrame.rot.coeffs();
	states.col(i).segment(1 + 7*(j) + 4, 3) = objectFrame.pos;
      }
      else
      {
	states.col(i).segment(1 + 7*(j) + 0, 3) = objectFrame.pos;
	states.col(i).segment(1 + 7*(j) + 3, 4) = objectFrame.rot.coeffs();
      }
    }
  }


  // WRITE OUTPUT FILE
  
  
  // now output the uncompressed state matrix
  igl::writeDMAT(output_d, states, ascii);
}
