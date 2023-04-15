#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>

template<typename ValueType>
class Trie
{
public:
    Trie();
    ~Trie() {resetHelper(root);}
    void reset();
    void insert(const std::string& key, const ValueType& value);
    std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;
    
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
private:
    struct tNode
    {
        char label;
        std::vector<ValueType> values;
        std::vector<tNode*> children;
    };
    tNode *root;
    
    void insertHelper(const std::string& key, tNode* cur, const ValueType& value);
    void resetHelper(tNode* cur);
    std::vector<ValueType> findHelperTrue(const std::string& key, tNode* cur) const;
    std::vector<ValueType> findHelperFalse(const std::string& key, tNode* cur) const;

};

template<typename ValueType>
inline
Trie<ValueType>::Trie()
{
    tNode *root_node = new tNode;   
    root_node->label = ' ';
    root = root_node;
}

template<typename ValueType>
void Trie<ValueType>::reset()
{
    if (root->children.size() == 0)
        return;
    resetHelper(root);
    tNode *new_node = new tNode;
    root = new_node;
}

template<typename ValueType>
void Trie<ValueType>::insert(const std::string& key, const ValueType& value)
{
    if (key.size() == 0)   
        return;
    insertHelper(key, root, value);
}

template<typename ValueType>
void Trie<ValueType>::insertHelper(const std::string& key, tNode* cur, const ValueType& value)
{
    
    if (key.size() == 0)
    {
        cur->values.push_back(value);
        return;
    }
    if (cur->children.size() == 0)
    {
        tNode *newnode = new tNode;
        newnode->label = key[0];
        cur->children.push_back(newnode);
        return insertHelper(key.substr(1), newnode, value);
    }
    else
    {
        for (typename std::vector<tNode*>::iterator it = cur->children.begin();
             it != cur->children.end(); ++it)
        {
            if (key[0] == (*it)->label)
                return insertHelper(key.substr(1), *it, value);
        }
        tNode *newnodeb = new tNode;
        newnodeb->label = key[0];
        cur->children.push_back(newnodeb);
        return insertHelper(key.substr(1), newnodeb, value);
    }
}

template<typename ValueType>
void Trie<ValueType>::resetHelper(tNode* cur)
{
    if (cur == nullptr)
        return;
    
    for (typename std::vector<tNode*>::iterator it = cur->children.begin();
         it != cur->children.end(); it++)
    {
        resetHelper(*it);
    }
    
    delete cur;
}


template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const
{
    std::vector<ValueType> find_values;
    
    if (exactMatchOnly == true)
        return findHelperTrue(key, root);
    else
    {
        for (typename std::vector<tNode*>::iterator it = root->children.begin();
             it != root->children.end(); ++it)
        {
            if (key[0] == (*it)->label)    
                return findHelperFalse(key.substr(1), *it);
        }
    }
    
    return find_values;   
    
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::findHelperTrue(const std::string& key, tNode* cur) const
{
    std::vector<ValueType> find_values;
    
    if (key.size() == 0)
    {
        for (typename std::vector<ValueType>::iterator it = cur->values.begin();
             it != cur->values.end(); ++it)
            find_values.push_back(*it);
        return find_values;
    }
    
    for (typename std::vector<tNode*>::iterator it = cur->children.begin();
         it != cur->children.end(); ++it)
    {
        if (key[0] == (*it)->label)
            return findHelperTrue(key.substr(1), *it);
    }
    
    return find_values;
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::findHelperFalse(const std::string& key, tNode* cur) const
{
    std::vector<ValueType> find_values;
    
    if (key.size() == 0)
    {
        std::vector<ValueType> new_find_values;
        for (typename std::vector<ValueType>::iterator it = cur->values.begin();
             it != cur->values.end(); ++it)
            new_find_values.push_back(*it);
        return new_find_values;
    }
    
    for (typename std::vector<tNode*>::iterator it = cur->children.begin();
         it != cur->children.end(); ++it)
    {
        if (key[0] == (*it)->label)
        {
            std::vector<ValueType> f_find_values = findHelperFalse(key.substr(1), *it);
            find_values.insert(find_values.begin(), f_find_values.begin(), f_find_values.end());
        }
        
        if (key[0] != (*it)->label)   
        {
            std::vector<ValueType> k_find_values = findHelperTrue(key.substr(1), *it);
            find_values.insert(find_values.begin(), k_find_values.begin(),
                               k_find_values.end());
        }
    }
    
    return find_values;
    
}

#endif  
