#ifndef IO
#define IO

#include "mesh.hpp"
#include "vtk_writer.hpp"

#include <string>
#include <map>
#include <memory>
#include <stdexcept>

class MeshWriter {
public:
    virtual ~MeshWriter() = default;
    virtual void write(const Mesh& mesh, const std::string& path) const = 0;
};

class MeshReader {
public:
    virtual ~MeshReader() = default;
    virtual Mesh read(const std::string& path) const = 0;
};

class VtkWriter : public MeshWriter {
public:
    void write(const Mesh& mesh, const std::string& path) const override;
};

class NastranWriter : public MeshWriter {
public:
    void write(const Mesh& mesh, const std::string& path) const override;
};

std::unique_ptr<MeshWriter> makeWriter(const std::string& format);

class GmshReader : public MeshReader {
public:
    Mesh read(const std::string& path) const override;
};

std::unique_ptr<MeshReader> makeReader(const std::string& format);

std::map<std::string, std::string> readConfig(const std::string& path);

#endif