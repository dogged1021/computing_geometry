/**
 * @file GeoComp.h
 * @author Zhang CHEN (zhang.chen@agile-robots.com)
 * @brief 空间几何计算
 * @version 0.1
 * @date 2023-07-10
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef GEOMETRIC_CALCULATION_H
#define GEOMETRIC_CALCULATION_H

#include <iostream>
#include <math.h>
#include <algorithm>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>

namespace ar
{
namespace geometricCalculation
{
    const double EPS = 0.000001;

    /**
     * @brief 平面方程 ax+by+cz+d=0
     *
     * @param rectangle 三个不共面点
     * @return true
     */
    bool eqPlanar(const std::vector<Eigen::Vector3d> &rectangle, double &a, double &b, double &c, double &d);

    /**
     * @brief 点point到面ax+by+cz+d=0的距离
     *
     * @param distance 距离
     * @return true
     */
    bool disP2Planar(const Eigen::Vector3d &point, const double &a, const double &b, const double &c, const double &d, double &distance);

    /**
     * @brief 计算长方形面积
     *
     * @param rectangles 三个顶点数组
     * @param square 面积
     * @return true
     */
    bool areaRectangle(const std::vector<Eigen::Vector3d> &rectangles, double square);

    /**
     * @brief 给定一条直线上两点point1,point2，直线外一点point0, 计算点在直线的投影V_point
     * @param point0
     * @param point1
     * @param point2
     * @param V_point
     * @return true
     */
    bool getVpoint(const Eigen::Vector3d &point0, const Eigen::Vector3d &point1, const Eigen::Vector3d &point2, Eigen::Vector3d V_point);

    /**
     * @brief 判断空间中两个长方形是否有重叠
     *
     * @param rectangle0 第一个长方形的三个顶点
     * @param rectangle1 第二个长方形的三个顶点
     * @return true
     */
    bool ifOverlap(const std::vector<Eigen::Vector3d> rectangle0, const std::vector<Eigen::Vector3d> rectangle1);

    //============================================================================================================
    bool eqPlanar(const std::vector<Eigen::Vector3d> &rectangles, double &a, double &b, double &c, double &d)
    {
        Eigen::Vector3d vec0 = rectangles[0] - rectangles[1];
        Eigen::Vector3d vec1 = rectangles[2] - rectangles[1];

        if ((fabs(vec0.dot(vec1)) - vec0.norm() * vec1.norm()) < EPS)
        {
            std::cout << "三点共线" << std::endl;
            return false;
        }
        else
        {
            a = rectangles[0].y() * (rectangles[1].z() - rectangles[2].z()) + rectangles[1].y() * (rectangles[2].z() - rectangles[0].z()) + rectangles[2].y() * (rectangles[0].z() - rectangles[1].z());
            b = rectangles[0].z() * (rectangles[1].x() - rectangles[2].x()) + rectangles[1].z() * (rectangles[2].x() - rectangles[0].x()) + rectangles[2].z() * (rectangles[0].x() - rectangles[1].x());
            c = rectangles[0].x() * (rectangles[1].y() - rectangles[2].y()) + rectangles[1].x() * (rectangles[2].y() - rectangles[0].y()) + rectangles[2].x() * (rectangles[0].y() - rectangles[1].y());
            d = -rectangles[0].x() * (rectangles[1].y() * rectangles[2].z() - rectangles[2].y() * rectangles[1].z()) - rectangles[1].x() * (rectangles[2].y() * rectangles[0].z() - rectangles[0].y() * rectangles[2].z()) - rectangles[2].x() * (rectangles[0].y() * rectangles[1].z() - rectangles[1].y() * rectangles[0].z());
        }
        return true;
    }

    bool disP2Planar(const Eigen::Vector3d &point, const double &a, const double &b, const double &c, const double &d, double &distance)
    {
        distance = fabs(a * point.x() + b * point.y() + c * point.z() + d) / sqrt(a * a + b * b + c * c);
        // std::cout << "点到面的距离为 " << distance << std::endl;
        return true;
    }

    bool areaRectangle(const std::vector<Eigen::Vector3d> &rectangles, double square)
    {
        Eigen::Vector3d l0, l1, l2;

        l0 = rectangles[0] - rectangles[1];
        l1 = rectangles[0] - rectangles[2];
        l2 = rectangles[1] - rectangles[2];
        if (fabs(l0.transpose() * l1) < EPS)
        {
            square = l0.norm() * l1.norm();
        }
        else if (fabs(l0.transpose() * l2) < EPS)
        {
            square = l0.norm() * l2.norm();
        }
        else
        {
            square = l1.norm() * l2.norm();
        }

        return true;
    }

    bool getVpoint(const Eigen::Vector3d &point0, const Eigen::Vector3d &point1, const Eigen::Vector3d &point2, Eigen::Vector3d V_point)
    {
        if (fabs(fabs((point0 - point1).dot(point2 - point1) / ((point0 - point1).norm() * (point2 - point1).norm())) - 1) < EPS)
        {
            V_point = point0;
        }
        else
        {
            std::vector<Eigen::Vector3d> rectangle0, rectangle1;
            rectangle0.push_back(point0);
            rectangle0.push_back(point1);
            rectangle0.push_back(point2);
            double a0, b0, c0, d0, a1, b1, c1, d1;
            eqPlanar(rectangle0, a0, b0, c0, d0);
            Eigen::Vector3d nvector2plan0(a0, b0, c0);
            Eigen::Vector3d point3 = point0 + nvector2plan0; // 新定义一个不在给点三点确定的平面上，以构造两个平面，其交线为p1p2
            rectangle1 = rectangle0;
            rectangle1[0] = point3;
            eqPlanar(rectangle1, a1, b1, c1, d1);

            // 求解方程组
            Eigen::Matrix3d A;
            Eigen::Vector3d b;

            A << a0, b0, c0,
                a1, b1, c1,
                (point1.x() - point2.x()), (point1.y() - point2.y()), (point1.z() - point2.z());
            b << -d0, -d1, point0.x() * (point1.x() - point2.x()) + point0.y() * (point1.y() - point2.y()) + point0.z() * (point1.z() - point2.z());
            V_point = A.lu().solve(b);
        }

        return true;
    }

    bool ifOverlap(const std::vector<Eigen::Vector3d> rectangle0, const std::vector<Eigen::Vector3d> rectangle1)
    {
        Eigen::Vector3d l0, l1, l2;
        Eigen::Vector3d point_center;
        Eigen::Vector3d point0[4], point1[4];

        // 按顺序确定顶点坐标
        // 计算第一个长方形的特征
        l0 = rectangle0[0] - rectangle0[1];
        l1 = rectangle0[0] - rectangle0[2];
        l2 = rectangle0[1] - rectangle0[2];

        if (fabs(l0.transpose() * l1) < EPS) // 第一个长方形的中心
        {
            point_center = (rectangle0[1] + rectangle0[2]) / 2; // 第一个长方形的中心
            point0[0] = rectangle0[1];
            point0[1] = rectangle0[0];
            point0[2] = rectangle0[2];
            point0[3] = 2 * point_center - rectangle0[0]; // 按包围顺序重新命名长方形四个顶点
        }
        else if (fabs(l0.transpose() * l2) < EPS)
        {
            point_center = (rectangle0[0] + rectangle0[2]) / 2;
            point0[3] = 2 * point_center - rectangle0[1];
            point0[0] = rectangle0[0];
            point0[1] = rectangle0[1];
            point0[2] = rectangle0[2];
        }
        else
        {
            point_center = (rectangle0[0] + rectangle0[1]) / 2;
            point0[3] = 2 * point_center - rectangle0[2];
            point0[0] = rectangle0[0];
            point0[1] = rectangle0[2];
            point0[2] = rectangle0[1];
        }

        // 计算第二个长方形的特征
        l0 = rectangle1[0] - rectangle1[1];
        l1 = rectangle1[0] - rectangle1[2];
        l2 = rectangle1[1] - rectangle1[2];

        if (fabs(l0.transpose() * l1) < EPS) // 第二个长方形的中心
        {
            point_center = (rectangle1[1] + rectangle1[2]) / 2; // 第二个长方形的中心
            point1[0] = rectangle1[1];
            point1[1] = rectangle1[0];
            point1[2] = rectangle1[2];
            point1[3] = 2 * point_center - rectangle1[0]; // 按包围顺序重新命名长方形四个顶点
        }
        else if (fabs(l0.transpose() * l2) < EPS)
        {
            point_center = (rectangle1[0] + rectangle1[2]) / 2;
            point1[3] = 2 * point_center - rectangle1[1];
            point1[0] = rectangle1[0];
            point1[1] = rectangle1[1];
            point1[2] = rectangle1[2];
        }
        else
        {
            point_center = (rectangle1[0] + rectangle1[1]) / 2;
            point1[3] = 2 * point_center - rectangle1[2];
            point1[0] = rectangle1[0];
            point1[1] = rectangle1[2];
            point1[2] = rectangle1[1];
        }

        // 投影
        Eigen::Vector3d ProjPoint[4];
        int flag[4] = {0};
        double MaxLength, SumLength;
        // 向第一个长方形第一条边上投影
        for (int j = 0; j < 4; j++)
        {
            getVpoint(point1[j], point0[0], point0[1], ProjPoint[j]);
        }

        MaxLength = std::max({(point0[0] - ProjPoint[0]).norm(), (point0[0] - ProjPoint[1]).norm(), (point0[0] - ProjPoint[2]).norm(), (point0[0] - ProjPoint[3]).norm(),
                              (point0[1] - ProjPoint[0]).norm(), (point0[1] - ProjPoint[1]).norm(), (point0[1] - ProjPoint[2]).norm(), (point0[1] - ProjPoint[3]).norm()});

        SumLength = (point0[0] - point0[1]).norm() +
                    std::max({(ProjPoint[0] - ProjPoint[1]).norm(), (ProjPoint[0] - ProjPoint[2]).norm(), (ProjPoint[0] - ProjPoint[3]).norm(),
                              (ProjPoint[1] - ProjPoint[2]).norm(), (ProjPoint[1] - ProjPoint[3]).norm(), (ProjPoint[2] - ProjPoint[3]).norm()});

        if ((MaxLength - SumLength) < EPS)
        {
            flag[0] = 1; // 此投影方向重叠
        }

        // 向第一个长方形第二条边上投影
        for (int j = 0; j < 4; j++)
        {
            getVpoint(point1[j], point0[1], point0[2], ProjPoint[j]);
        }

        MaxLength = std::max({(point0[2] - ProjPoint[0]).norm(), (point0[2] - ProjPoint[1]).norm(), (point0[2] - ProjPoint[2]).norm(), (point0[2] - ProjPoint[3]).norm(),
                              (point0[1] - ProjPoint[0]).norm(), (point0[1] - ProjPoint[1]).norm(), (point0[1] - ProjPoint[2]).norm(), (point0[1] - ProjPoint[3]).norm()});
        SumLength = (point0[2] - point0[1]).norm() +
                    std::max({(ProjPoint[0] - ProjPoint[1]).norm(), (ProjPoint[0] - ProjPoint[2]).norm(), (ProjPoint[0] - ProjPoint[3]).norm(),
                              (ProjPoint[1] - ProjPoint[2]).norm(), (ProjPoint[1] - ProjPoint[3]).norm(), (ProjPoint[2] - ProjPoint[3]).norm()});

        if ((MaxLength - SumLength) < EPS)
        {
            flag[1] = 1; // 此投影方向重叠
        }
        // 向第二个长方形第一条边上投影
        for (int j = 0; j < 4; j++)
        {
            getVpoint(point0[j], point1[0], point1[1], ProjPoint[j]);
        }
        MaxLength = std::max({(point1[0] - ProjPoint[0]).norm(), (point1[0] - ProjPoint[1]).norm(), (point1[0] - ProjPoint[2]).norm(), (point1[0] - ProjPoint[3]).norm(),
                              (point1[1] - ProjPoint[0]).norm(), (point1[1] - ProjPoint[1]).norm(), (point1[1] - ProjPoint[2]).norm(), (point1[1] - ProjPoint[3]).norm()});
        SumLength = (point1[0] - point1[1]).norm() +
                    std::max({(ProjPoint[0] - ProjPoint[1]).norm(), (ProjPoint[0] - ProjPoint[2]).norm(), (ProjPoint[0] - ProjPoint[3]).norm(),
                              (ProjPoint[1] - ProjPoint[2]).norm(), (ProjPoint[1] - ProjPoint[3]).norm(), (ProjPoint[2] - ProjPoint[3]).norm()});
        // std::cout << "MaxLength = " << MaxLength << std::endl << "SumLength = " << SumLength << std::endl;
        if ((MaxLength - SumLength) < EPS)
        {
            flag[2] = 1; // 此投影方向重叠
        }

        for (int j = 0; j < 4; j++)
        {
            getVpoint(point0[j], point1[1], point1[2], ProjPoint[j]);
        }
        MaxLength = std::max({(point1[2] - ProjPoint[0]).norm(), (point1[2] - ProjPoint[1]).norm(), (point1[2] - ProjPoint[2]).norm(), (point1[2] - ProjPoint[3]).norm(),
                              (point1[1] - ProjPoint[0]).norm(), (point1[1] - ProjPoint[1]).norm(), (point1[1] - ProjPoint[2]).norm(), (point1[1] - ProjPoint[3]).norm()});
        SumLength = (point1[2] - point1[1]).norm() +
                    std::max({(ProjPoint[0] - ProjPoint[1]).norm(), (ProjPoint[0] - ProjPoint[2]).norm(), (ProjPoint[0] - ProjPoint[3]).norm(),
                              (ProjPoint[1] - ProjPoint[2]).norm(), (ProjPoint[1] - ProjPoint[3]).norm(), (ProjPoint[2] - ProjPoint[3]).norm()});
        // std::cout << "MaxLength = " << MaxLength << std::endl << "SumLength = " << SumLength << std::endl;
        if ((MaxLength - SumLength) < EPS)
        {
            flag[3] = 1; // 此投影方向重叠
        }

        if (flag[0] == 1 && flag[1] == 1 && flag[2] == 1 && flag[3] == 1)
        {
            return true; // 每一个投影方向都重叠，则重叠
        }
        else
        {
            return false; // 存在不重叠的投影方向，则不重叠
        }
    }
} // namespace geometricCalculation
} // namespace ar
#endif