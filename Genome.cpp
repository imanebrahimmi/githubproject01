#include "Animall.h"
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <istream>

using namespace std;

class GenomeImpl
{
public:
    GenomeImpl(const string& nm, const string& sequence);
    static bool load(istream& genomeSource, vector<Genome>& genomes);
    int length() const;
    string name() const;
    bool extract(int position, int length, string& fragment) const;
private:
    string m_name;
    string m_sequence;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
{
    m_name = nm;
    m_sequence = sequence;
}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes)
{
    string s;
    string name;
    string dna;
    string dna_upper;
    bool isfirst = true;
    while (getline(genomeSource, s))
    {
        if (s.at(0) == '>')
        {
            if (!isfirst) {
                for (int i = 0; i < dna.size(); i++)
                {
                    if (dna[i] == 'a' || dna[i] == 'A')
                        dna_upper += "A";
                    if (dna[i] == 'c' || dna[i] == 'C')
                        dna_upper += "C";
                    if (dna[i] == 't' || dna[i] == 'T')
                        dna_upper += "T";
                    if (dna[i] == 'g' || dna[i] == 'G')
                        dna_upper += "G";
                    if (dna[i] == 'n' || dna[i] == 'N')
                        dna_upper += "N";
                }
                Genome new_genome(name, dna_upper);
                genomes.push_back(new_genome);
            }
            else isfirst = false;
            
            name = s.substr(1);   
            dna = "";
            dna_upper = "";
        }
        else
        {
            dna += s;   
        }
    }
    
    for (int i = 0; i < dna.size(); i++)
    {
        if (dna[i] == 'a' || dna[i] == 'A')
            dna_upper += "A";
        if (dna[i] == 'c' || dna[i] == 'C')
            dna_upper += "C";
        if (dna[i] == 't' || dna[i] == 'T')
            dna_upper += "T";
        if (dna[i] == 'g' || dna[i] == 'G')
            dna_upper += "G";
        if (dna[i] == 'n' || dna[i] == 'N')
            dna_upper += "N";
    }
    Genome new_genome(name, dna_upper);  
    genomes.push_back(new_genome);
    
    return true;
}

int GenomeImpl::length() const
{
    return (int)m_sequence.size();
}

string GenomeImpl::name() const
{
    return m_name;
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
    if (position > m_sequence.size())
        return false;
    if (position + length > m_sequence.size())
        return false;
    
    for (int i = position; i < (position + length); i++)
    {
        fragment += m_sequence[i];
    }
    return true;
}
 

Genome::Genome(const string& nm, const string& sequence)
{
    m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
    delete m_impl;
}

Genome::Genome(const Genome& other)
{
    m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
    GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
    delete m_impl;
    m_impl = newImpl;
    return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes)
{
    return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
    return m_impl->length();
}

string Genome::name() const
{
    return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
    return m_impl->extract(position, length, fragment);
}
