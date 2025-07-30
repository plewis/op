#pragma once

namespace op
{

    class Tree;
    class TreeManip;
    class Likelihood;
    //class Updater;

    class Node {
            friend class Tree;
            friend class TreeManip;
            friend class Likelihood;
            //friend class Updater;

        public:
                                        Node();
                                        ~Node();

                    Node *              getParent()     {return _parent;}
                    Node *              getLeftChild()  {return _left_child;}
                    Node *              getRightSib()   {return _right_sib;}
                    int                 getNumber()     {return _number;}
                    string              getName()       {return _name;}
                    Split               getSplit()      {return _split;}

                    double              getEdgeLength() {return _edge_length;}
                    void                setEdgeLength(double v);

            static const double _smallest_edge_length;

            typedef vector<Node>    Vector;
            typedef vector<Node *>  PtrVector;

        private:

            void                clear();

            Node *              _left_child;
            Node *              _right_sib;
            Node *              _parent;
            int                 _number;
            string              _name;
            double              _edge_length;
            Split               _split;
    };

    inline Node::Node() {
        //cout << "Creating Node object" << endl;
        clear();
    }

    inline Node::~Node() {
        //cout << "Destroying Node object" << endl;
    }

    inline void Node::clear() {
        _left_child = nullptr;
        _right_sib = nullptr;
        _parent = nullptr;
        _number = 0;
        _name = "";
        _edge_length = _smallest_edge_length;
        _split.clear();
    }

    inline void Node::setEdgeLength(double v) {
        _edge_length = (v < _smallest_edge_length ? _smallest_edge_length : v);
    }

}

