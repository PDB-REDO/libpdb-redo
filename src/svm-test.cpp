/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "pdb-redo.hpp"

#include <fstream>

#include <zeep/xml/serialize.hpp>
#include <zeep/xml/document.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "svm++.h"

namespace ba = boost::algorithm;

//struct svm_class
//{
//	int8_t label;
//	size_t nr_sv;
//	
//	template<class Archive>
//	void serialize(Archive& ar, const unsigned int version)
//	{
//		ar & zeep::xml::make_attribute_nvp("label", label)
//		   & zeep::xml::make_attribute_nvp("nr_sv", nr_sv);
//	}
//};
//
//struct svm_sv
//{
//	int index;
//	double value;
//
//	template<class Archive>
//	void serialize(Archive& ar, const unsigned int version)
//	{
//		ar & zeep::xml::make_attribute_nvp("index", index)
//		   & zeep::xml::make_element_nvp(".", value);
//	}
//};
//
//struct svm_node
//{
//	std::vector<double> sv_coef;
//	std::vector<svm_sv> sv;
//
//	template<class Archive>
//	void serialize(Archive& ar, const unsigned int version)
//	{
//		ar & zeep::xml::make_element_nvp("sv_coef", sv_coef)
//		   & zeep::xml::make_element_nvp("sv", sv);
//	}
//};
//
//struct svm_config
//{
//	std::string					svm_type;
//	std::string					kernel_type;
//	boost::optional<double>	gamma;
//	std::vector<double>			rho;
//	std::vector<svm_class>		classes;
//	std::vector<svm_node>		sv;
//
//	template<class Archive>
//	void serialize(Archive& ar, const unsigned int version)
//	{
//		ar & zeep::xml::make_element_nvp("svm-type", svm_type)
//		   & zeep::xml::make_element_nvp("kernel-type", kernel_type)
//		   & zeep::xml::make_element_nvp("gamma", gamma)
//		   & zeep::xml::make_element_nvp("rho", rho)
//		   & zeep::xml::make_element_nvp("class", classes)
//		   & zeep::xml::make_element_nvp("svm-node", sv);
//	}
//};

void predict_file(svm::ModelBase* model, const char* filename)
{
//	typedef svm::ModelBase	SVMModel;
	typedef svm::Vector		SVMVector;

	std::ifstream test(filename);
	std::string line;

	while (std::getline(test, line))
	{
		std::list<std::string> f;
		ba::split(f, line, ba::is_any_of(" "));
		f.pop_front();
		
		SVMVector v;
		for (std::string& fs: f)
		{
			auto c = fs.find(':');
			if (c == std::string::npos)
				continue;
			
			int index = boost::lexical_cast<int>(fs.substr(0, c));
			float value = boost::lexical_cast<float>(fs.substr(c + 1));
			
			v[index] = value;
		}

		char predictedClass = model->Predict(v);
	
		auto prob = model->PredictWithProbability(v);

		auto voted = distance(prob.begin(),
			max_element(prob.begin(), prob.end()));
		
		std::cout << predictedClass << '\t'
			 << (voted ? '0' : '1') << '\t'
			 << "probabilities: ";
		copy(prob.begin(), prob.end(), std::ostream_iterator<double>(std::cout, ", "));
		std::cout << std::endl;
	}
}

//int centrifuge_predict(int argc, char* argv[])
//{
//	ifstream file("1cbs-0.8.svm");
//	
////	svm_config config;
//	
//	zeep::xml::document doc(file);
////	doc.deserialize("svm-config", config);
//
//	auto config = doc.find_first("//svm-config");
//	if (not config)
//		throw std::runtime_error("invalid svm file");
//	
//	auto model = svm::ModelBase::Create(*config);
//
//	predict_file(model, "1ctn-0.8-scaled.txt");
//
//	return 0;
//}

int centrifuge_test(int argc, char* argv[])
{
	cif::VERBOSE = 1;
	
	std::ifstream data("1cbs-0.8.txt.scale");
	if (not data.is_open())
		throw std::runtime_error("no such file");

	typedef svm::SVM_C_SVC_RBF SVM;
	typedef SVM::param_type SVMParams;
	typedef svm::Matrix SVMMatrix;

	SVMParams params;
	params.gamma = 0.5;
	params.C = 8;
	params.probability = true;

	SVM svm(params);

	std::vector<int8_t> labels;
	SVMMatrix m;
	size_t row = 0;

	for (;;)
	{
		std::string line;
		std::getline(data, line);
		
		if (line.empty())
		{
			if (data.eof())
				break;
			continue;
		}
		
		std::vector<std::string> f;
		ba::split(f, line, ba::is_any_of(" "));
		
		if (f.size() < 1)
			throw std::runtime_error("invalid data");
		
		int8_t c = boost::lexical_cast<int8_t>(f.front());
		f.erase(f.begin(), f.begin() + 1);
		
		labels.push_back(c);
		
		for (auto& fs: f)
		{
			auto c = fs.find(':');
			if (c == std::string::npos)
				continue;
			
			int index = boost::lexical_cast<int>(fs.substr(0, c));
			float value = boost::lexical_cast<float>(fs.substr(c + 1));

			m[row][index - 1] = value;
		}
		
		++row;
	}

	auto model = svm.Train(labels, m);
	
//	model->Print();

	auto cv = svm.CrossValidation(labels, m, 5);
	
	size_t correct = 0;
	for (size_t i = 0; i < labels.size(); ++i)
		if (cv[i] == labels[i])
			++correct;
	
	std::cout << "Cross validation Accuracy = " << (100.0 * correct / labels.size()) << '%' << std::endl;

	predict_file(model, "1ctn-0.8-scaled.txt");
	
	zeep::xml::document doc;
	doc.root()->append(model->GetConfig());
	
	std::ofstream model_file("test.model");
	if (model_file.is_open())
		model_file << doc;

	return 0;
}
