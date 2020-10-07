//
// Created by Sultan Sinem Eren on 18.07.20.
//

#ifndef CODINGTASK_INTERSECT_H
#define CODINGTASK_INTERSECT_H

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <math.h>


#define PRINT(x) std::cout << x << std::endl;

class Intersect {
public:
    static bool overlap(std::string polygonA, std::string polygonB);

private:
};

struct Point{
    double x;
    double y;
};

/// parse the vertices from the string
/// @param <polygon> polygon vertices as a string
/// @return vertices of the polygon
std::vector<Point> getVertices(const std::string &polygon)
{

    std::vector<Point> vertices;
    vertices.reserve(polygon.size());

    std::stringstream s_stream(polygon);
    std::string::size_type sz;
    Point point;

    while(s_stream.good())
    {
        std::string substr;
        getline(s_stream, substr, ',');

        point.x = std::stod(substr,&sz);
        point.y = std::stod(substr.substr(sz));

        vertices.push_back(point);
    }


    return vertices;
}

/// find the edges of the polygon
/// @param <vertices> vertices of the polygon
/// @return edges of the polygon
std::vector<Point> getEdges(const std::vector<Point> &vertices)
{
    std::vector<Point> edges;
    edges.reserve(vertices.size());
    Point point;

    for(size_t i = 0; i < vertices.size()-1; ++i)
    {
        point.x = vertices.at(i+1).x-vertices.at(i).x;
        point.y = vertices.at(i+1).y - vertices.at(i).y;

        edges.push_back(point);
    }



    return edges;


}

/// find the normals of each edge of the polygon
/// @param <edges> edges of the polygon
/// @return normals of each edge of the polygon
std::vector<Point> getPerpendicularVect(const std::vector<Point> &edges)
{

    std::vector<Point> perpendicularVects;
    perpendicularVects.reserve(edges.size());

    Point point;

    for(auto &edge:edges)
    {
        point.x = -edge.y;
        point.y = edge.x;

        perpendicularVects.push_back(point);
    }


    return perpendicularVects;
}

/// order the vertices in clockwise
/// @param <vertices> vertices of the polygon
/// @return clockwise ordered vertices
std::vector<Point> ClokwiseOrder(std::vector<Point> &vertices)
{

    double center_x = 0.0;
    double center_y = 0.0;


    std::vector<std::pair<double,double>> thetaAndIndex;
    thetaAndIndex.reserve(vertices.size());
    std::vector<Point> new_vertices(vertices.size());
    double theta;


    for(auto &vertex:vertices)
    {
        center_x+=vertex.x;
        center_y+=vertex.y;
    }

    center_x/=vertices.size();
    center_y/=vertices.size();

    for (size_t i = 0; i < vertices.size()-1; ++i)
    {
        theta = ((atan2(vertices.at(i).x - center_x, vertices.at(i).y - center_y) * 180.0 / M_PI )+ 360.0);

        thetaAndIndex.push_back(std::make_pair(theta,i));
    }


    std::sort(thetaAndIndex.begin(),thetaAndIndex.end());

    for(size_t i = 0; i < vertices.size()-1; ++i)
    {
        new_vertices.at(i) = vertices.at(thetaAndIndex.at(i).second);

    }

    new_vertices.at(vertices.size()-1) = vertices.at(thetaAndIndex.at(0).second);


    return new_vertices;

}

/// convexixt check
/// only works if polygon is not self intersecting
/// @param <vertices> vertices of the polygon
/// @return if the polygon is convex it returns true
bool isConvex(std::vector<Point> vertices)
{
    // if the polygon is triangle
    if(vertices.size()-1 > 3)
    {
        return true;
    }

    double cross, eX1, eX2, eY1, eY2;
    int num = vertices.size();
    bool signOfZ = false;

    for( size_t i = 0; i < num-1; ++i)
    {

        eX1 = vertices.at((i+2) % num).x - vertices.at((i+1)% num).x;
        eY1 = vertices.at((i+2) % num).y - vertices.at((i+1)% num).y;

        eX2 = vertices.at(i).x - vertices.at((i+1) % num).x;
        eY2 = vertices.at(i).y - vertices.at((i+1) % num).y;

        cross = eX1 * eY2 - eY1 * eX2;

        if(i == 0)
        {
            signOfZ = (cross > 0.0);
        }
        else if( signOfZ != (cross>0.0))
        {
            return false;
        }


    }

    return true;

}

/// This function checks whether the polygon is closed or not
/// It throw a logic_error if the polygon is not closed
/// @param <vertices> vertices of the polygon
void isClosed(std::vector<Point> vertices)
{

    bool xCheck;
    bool yCheck;
    unsigned int n = vertices.size()-1;

    xCheck = (vertices.at(0).x == vertices.at(n).x );
    yCheck = (vertices.at(0).y == vertices.at(n).y);


    if( ((!xCheck)  || (!yCheck)) )
    {
        throw std::logic_error(" THE POLYGON IS NOT CLOSED!");
    }


}

/// determine whether two polygons overlap or not
/// based on seperating axis theorem so it only works for convex polygons
/// @param <polygonA> first polygon
/// @param <polygonB> second polygon
/// @return true if polygons overlap and false if they don't overlap
bool Intersect::overlap(std::string polygonA, std::string polygonB) {

    // your implementation goes here

    /*
     * Second task answer:
     * No. It does not work for non-convex polygons because I used the separating axis theorem.
     * Non-convex polygons can be represented by combination of convex polygons.
     * If the convex decomposition is implemented then we can devide
     * the non-convex polygons into convex polygons and run the  SAT.
     */

    // vertices of the polygons
    std::vector<Point> vertices_A = getVertices(polygonA);
    std::vector<Point> vertices_B = getVertices(polygonB);

    isClosed(vertices_B);
    isClosed(vertices_A);


    // to make sure the clockwise ordering
    vertices_A = ClokwiseOrder(vertices_A);
    vertices_B = ClokwiseOrder(vertices_B);



    // edges of the polygons
    std::vector<Point> edges_A = getEdges(vertices_A);
    std::vector<Point> edges_B = getEdges(vertices_B);



    // get perpendicular lines for all edges (normals)
    std::vector<Point> perpendicularVects_A = getPerpendicularVect(edges_A);
    std::vector<Point> perpendicularVects_B = getPerpendicularVect(edges_B);


    // check for the vertices
    bool xCheck, yCheck;
    for(size_t i = 0; i < vertices_A.size(); ++i)
    {
        for(size_t j = 0; j < vertices_B.size(); ++j)
        {
            xCheck = (vertices_A.at(i).x == vertices_B.at(j).x);
            yCheck = (vertices_A.at(i).y == vertices_B.at(j).y);

            if(xCheck && yCheck)
            {
                return true;
            }
        }
    }

    double A_min, A_max, B_min, B_max;
    double projection;

    //normal of polygon A
    for( size_t i = 0; i < perpendicularVects_A.size(); ++i)
    {


        // setting initial min and max
        A_min =  perpendicularVects_A[i].x * vertices_A[0].x + perpendicularVects_A[i].y * vertices_A[0].y;
        A_max = A_min;

        for( size_t j = 1; j < vertices_A.size(); ++j)
        {
            projection = perpendicularVects_A[i].x * vertices_A[j].x + perpendicularVects_A[i].y * vertices_A[j].y;

            if( projection < A_min)
            {
                A_min = projection;
            }
            else if( projection > A_min)
            {
                A_max = projection;
            }
        }


        // setting initial min and max
        B_min =  perpendicularVects_A[i].x * vertices_B[0].x + perpendicularVects_A[i].y * vertices_B[0].y;
        B_max = B_min;

        for( size_t j = 1; j < vertices_B.size(); ++j)
        {
            projection = perpendicularVects_A[i].x * vertices_B[j].x + perpendicularVects_A[i].y * vertices_B[j].y;

            if( projection < B_min)
            {
                B_min = projection;
            }
            else if( projection > B_min)
            {
                B_max = projection;
            }
        }

        if (!(A_min <= B_max && A_max >= B_min))
        {
            return false;
        }

    }




    //normal of polygon B
    for( size_t i = 0; i < perpendicularVects_B.size(); ++i)
    {



        // setting initial min and max for polygon A
        A_min =  perpendicularVects_B[i].x * vertices_A[0].x + perpendicularVects_B[i].y * vertices_A[0].y;
        A_max = A_min;

        for( size_t j = 1; j < vertices_A.size(); ++j)
        {
            projection = perpendicularVects_B[i].x * vertices_A[j].x + perpendicularVects_B[i].y * vertices_A[j].y;

            if( projection < A_min)
            {
                A_min = projection;
            }
            else if( projection > A_min)
            {
                A_max = projection;
            }
        }


        // setting initial min and max for polygon B
        B_min =  perpendicularVects_A[i].x * vertices_B[0].x + perpendicularVects_A[i].y * vertices_B[0].y;
        B_max = B_min;

        for( size_t j = 1; j < vertices_B.size(); ++j)
        {
             projection = perpendicularVects_B[i].x * vertices_B[j].x + perpendicularVects_B[i].y * vertices_B[j].y;

            if( projection < B_min)
            {
                B_min = projection;
            }
            else if( projection > B_min)
            {
                B_max = projection;
            }
        }

        if (!(A_min < B_max && A_max > B_min))
        {
            return false;
        }

    }


    return true;

}



#endif //CODINGTASK_INTERSECT_H
