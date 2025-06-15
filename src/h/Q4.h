/*****************************************************************************/
/*  STAP++ : CQ4 4节点平面单元（支持平面应力/应变分析）                       */
/*  参考Abaqus实现，支持结构平面问题建模                                      */
/*                                                                           */
/*  Computational Dynamics Laboratory                                        */
/*  School of Aerospace Engineering, Tsinghua University                     */
/*                                                                           */
/*  Release 1.11, November 22, 2017                                          */
/*                                                                           */
/*  http://www.comdyn.cn/                                                   */
/*****************************************************************************/

#pragma once

#include "Element.h"

// CQ4平面四节点单元类
class CQ4 : public CElement {
public:
    // 构造函数
    CQ4();

    // 析构函数
    ~CQ4();

    // 读取单元数据
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

    // 输出单元数据
    virtual void Write(COutputter& output);

    // 计算单元刚度矩阵
    virtual void ElementStiffness(double* Matrix);

    // 计算单元应力
    virtual void ElementStress(double* stress, double* Displacement);

    // 获取第i个节点
    CNode* GetNode(unsigned int i) const;

private:
    // 可扩展的单元参数（如厚度、应力状态标志等）
};
