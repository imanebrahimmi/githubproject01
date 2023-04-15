#ifndef PROVIDED_INCLUDED
#define PROVIDED_INCLUDED

#include <string>
#include <vector>
#include <istream>

class GenomeImpl;

class Genome
{
public:
    Genome(const std::string& nm, const std::string& sequence);
    ~Genome();
    Genome(const Genome& other);
    Genome& operator=(const Genome& rhs);
    static bool load(std::istream& genomeSource, std::vector<Genome>& genomes);
    int length() const;
    std::string name() const;
    bool extract(int position, int length, std::string& fragment) const;
    
private:
    GenomeImpl* m_impl;
};

struct DNAMatch
{
    std::string genomeName;
    int length;
    int position;
};

struct GenomeMatch
{
    std::string genomeName;
    double percentMatch;
};

class animallImpl;

class animall
{
public:
    animall(int minSearchLength);
    ~animall();
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const std::string& fragment, int minimumLength, bool exactMatchOnly, std::vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, std::vector<GenomeMatch>& results) const;
    
    animall(const animall&) = delete;
    animall& operator=(const animall&) = delete;
    
private:
    animallImpl* m_impl;
};

#endif 
