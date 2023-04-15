#include "Animall.h"
#include "Cell.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

class animallImpl
{
public:
    animallImpl(int minSearchLength);
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
    int m_minSearchLength; 
    
    struct GenomeAndPos
    {
        string genome_name;
        int pos;
    };
    
    vector<Genome> m_genome;
    Trie<GenomeAndPos> trie;
    void whatPushBack(vector<DNAMatch>& matches, DNAMatch new_match) const;
    void findRelatedGenomesHelper(vector<GenomeMatch>& temp_results, GenomeMatch new_GenomeMatch) const;
};

animallImpl::animallImpl(int minSearchLength)
{
    m_minSearchLength = minSearchLength;
}

void animallImpl::addGenome(const Genome& genome)
{
    string whole_sequence;
    genome.extract(0, genome.length(), whole_sequence);
    
    m_genome.push_back(genome);
    
    if (minimumSearchLength() >= genome.length())
    {
        string dna_sequency;
        if (genome.extract(0, minimumSearchLength(), dna_sequency))
        {
            GenomeAndPos temp;
            temp.genome_name = genome.name();
            temp.pos = 0;
            trie.insert(dna_sequency, temp);
        }
    }
    else
    {
        for (int i = 0; i <= (genome.length() - minimumSearchLength()); i++)
        {
            string dna_sequency;
            if (genome.extract(i, minimumSearchLength(), dna_sequency))
            {
                GenomeAndPos temp;
                temp.genome_name = genome.name();
                temp.pos = i;
                trie.insert(dna_sequency, temp);
            }
        }
    }
}

int animallImpl::minimumSearchLength() const
{
    return m_minSearchLength;
}

void animallImpl::whatPushBack(vector<DNAMatch>& matches, DNAMatch new_match) const
{
    if (matches.size() == 0)
        matches.push_back(new_match);
    
    bool temp1 = false; 
    bool temp2 = false; 
    
    for (vector<DNAMatch>::iterator it = matches.begin(); it != matches.end(); ++it)
    {
        if ((*it).length == new_match.length)
        {
            if ((*it).genomeName == new_match.genomeName)
                temp1 = true;
        }
        if ((*it).genomeName == new_match.genomeName)
        {
            if ((*it).length < new_match.length)
            {
                (*it).position = new_match.position;
                (*it).length = new_match.length;
                temp2 = true;
            }
            else
                return;
        }
    }
    if (!temp1 && !temp2)
        matches.push_back(new_match);
}


bool animallImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    if (fragment.size() < minimumLength)
        return false;
    if (minimumLength < minimumSearchLength())
        return false;
    
    string prefi_sequence;
    bool extract_bool = true;
    int length = 0;
    int misMatches = 0;
    for (int i = 0; i < minimumSearchLength(); i++)
    {
        prefi_sequence += fragment[i];
    }
    
    if (fragment.size() == minimumSearchLength())
    {
        vector<GenomeAndPos> temp = trie.find(prefi_sequence, exactMatchOnly);
        for (vector<GenomeAndPos>::iterator it = temp.begin(); it != temp.end(); ++it)
        {
            DNAMatch new_match;
            new_match.genomeName = (*it).genome_name;
            new_match.position = (*it).pos;
            new_match.length = minimumLength;
            whatPushBack(matches, new_match);
        }
    }
    
    if (fragment.size() > minimumSearchLength())
    {
        vector<GenomeAndPos> temp = trie.find(prefi_sequence, exactMatchOnly);
        for (vector<GenomeAndPos>::iterator it = temp.begin(); it != temp.end(); ++it)
        {
            string match_sequence;
            for (int i = 0; i < m_genome.size(); i++)
            {
                if (m_genome[i].name() == (*it).genome_name)
                {
                    if(((*it).pos + fragment.size()) < m_genome[i].length())
                        extract_bool = m_genome[i].extract((*it).pos, fragment.size(), match_sequence);
                    else if(((*it).pos + fragment.size()) > m_genome[i].length())
                        extract_bool = m_genome[i].extract((*it).pos, m_genome[i].length()-(*it).pos, match_sequence);
                }
            }
            if (extract_bool == false) 
                return false;
            
            for (int i = 0; i < minimumLength; i++)  
            {
                if (exactMatchOnly == true)
                {
                    if (match_sequence[i] == fragment[i])
                        length++;
                    else
                    {
                        length++;
                        misMatches++;
                    }
                }
                else if (exactMatchOnly == false)
                {
                    if (match_sequence[i] == fragment[i])
                        length++;
                    else if (match_sequence[i] != fragment[i])
                    {
                        length++;
                        misMatches++;
                    }
                }
            }
            
            if (exactMatchOnly == true)
            {
                if (misMatches < 1)
                {
                    for (int j = minimumLength; j < match_sequence.size(); j++)   
                    {
                        if (exactMatchOnly == true)
                        {
                            if (match_sequence[j] == fragment[j])
                                length++;
                            else
                                break;
                        }
                        else if (exactMatchOnly == false)
                        {
                            if (match_sequence[j] == fragment[j])
                                length++;
                            else if (match_sequence[j] != fragment[j])
                                break;
                        }
                    }
                }
            }
            if (exactMatchOnly == false)
            {
                if (misMatches < 2)
                {
                    for (int j = minimumLength; j < match_sequence.size(); j++)  
                    {
                        if (exactMatchOnly == true)
                        {
                            if (match_sequence[j] == fragment[j])
                                length++;
                            else
                                break;
                        }
                        else if (exactMatchOnly == false)
                        {
                            if (match_sequence[j] == fragment[j])
                                length++;
                            else if (match_sequence[j] != fragment[j])
                            {
                                if (misMatches == 0)
                                    length++;
                                else if (misMatches == 1)
                                    break;
                            }
                        }
                    }
                }
            }
            
            DNAMatch new_match;
            new_match.genomeName = (*it).genome_name;
            new_match.position = (*it).pos;
            new_match.length = length;
            if (exactMatchOnly == true && misMatches < 1)
                whatPushBack(matches, new_match);
            
            if (exactMatchOnly == false && misMatches < 2)
                whatPushBack(matches, new_match);
            
            length = 0;  
            misMatches = 0;
            extract_bool = true;
        }
        
    }
    
    if (matches.size() != 0)
        return true;
    else
        return false;
}

void animallImpl::findRelatedGenomesHelper(vector<GenomeMatch>& temp_results, GenomeMatch new_GenomeMatch) const
{
    if (temp_results.size() == 0)
    {
        temp_results.push_back(new_GenomeMatch);
        return;
    }
    
    if (temp_results.size() != 0)
    {
        for (vector<GenomeMatch>::iterator it = temp_results.begin(); it != temp_results.end(); ++it)
        {
            if ((*it).genomeName == new_GenomeMatch.genomeName)
            {
                (*it).percentMatch += new_GenomeMatch.percentMatch;
                return;
            }
        }
        
        temp_results.push_back(new_GenomeMatch);
    }
}

bool animallImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    if (fragmentMatchLength < minimumSearchLength())
        return false;
    
    vector<GenomeMatch> temp_results;
    bool extract_success;
    int divideNum = query.length()/fragmentMatchLength;
    for (int i = 0; i < divideNum; i++)
    {
        string sub_sequence;
        
        extract_success = query.extract(i*fragmentMatchLength, fragmentMatchLength, sub_sequence);
        if (extract_success == true)
        {
            vector<DNAMatch> matches;
            if (findGenomesWithThisDNA(sub_sequence, fragmentMatchLength, exactMatchOnly, matches))
            {
                for (int j = 0; j < matches.size(); j++)
                {
                    GenomeMatch new_GenomeMatch;
                    new_GenomeMatch.genomeName = matches[j].genomeName;
                    new_GenomeMatch.percentMatch = 1.0 / divideNum * 100.0;
                    findRelatedGenomesHelper(temp_results, new_GenomeMatch);
                }
            }
        }
    }
    
    for (int u = 0; u < temp_results.size(); u++)
    {
        if(temp_results[u].percentMatch >= matchPercentThreshold/100)
            results.push_back(temp_results[u]);
    }
    
    if (results.size() != 0)
        return true;
    else
        return false;
}

animall::animall(int minSearchLength)
{
    m_impl = new animallImpl(minSearchLength);
}

animall::~animall()
{
    delete m_impl;
}

void animall::addGenome(const Genome& genome)
{
    m_impl->addGenome(genome);
}

int animall::minimumSearchLength() const
{
    return m_impl->minimumSearchLength();
}

bool animall::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool animall::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}
