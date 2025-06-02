/*****************************************************************************/
/*  STAP++ : S4R 4节点壳单元（减少积分，含hourglass控制）                    */
/*  参考Abaqus实现，支持与其他3D单元混用                                      */
/*****************************************************************************/

#pragma once

#include "Element.h"

// S4R壳单元类
class CS4R : public CElement {
public:
    CS4R();
    ~CS4R();

    // 读取单元数据
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList) override;
    // 输出单元数据
    virtual void Write(COutputter& output) override;
    // 计算单元刚度矩阵（含hourglass控制）
    virtual void ElementStiffness(double* Matrix) override;
    // 计算单元应力
    virtual void ElementStress(double* stress, double* Displacement) override;

private:
    // hourglass 控制参数等可在此扩展
};
