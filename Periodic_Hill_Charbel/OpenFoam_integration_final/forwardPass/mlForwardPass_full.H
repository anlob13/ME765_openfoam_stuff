#ifndef forwardPass_full_H
#define forwardPass_full_H

#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>

#include "scalar.H"
#include "fileName.H"
#include "error.H"
#include <memory>
#include <utility>

namespace Foam
{

class Layer
{
public:

    virtual ~Layer() = default;

    virtual std::vector<scalar> forward
    (
        const std::vector<scalar>& x
    ) const = 0;
};


class DenseLayer : public Layer
{
public:

    label inDim;
    label outDim;

    std::vector<scalar> W;
    std::vector<scalar> b;

    scalar (*activation)(scalar);

    DenseLayer
    (
        label inDim_,
        label outDim_,
        std::vector<scalar> W_,
        std::vector<scalar> b_,
        scalar (*act)(scalar)
    )
    :
        inDim(inDim_),
        outDim(outDim_),
        W(std::move(W_)),
        b(std::move(b_)),
        activation(act)
    {}

    std::vector<scalar> forward
    (
        const std::vector<scalar>& x
    ) const override;
};


class NormalizationLayer : public Layer
{
public:

    std::vector<scalar> mean;
    std::vector<scalar> variance;

    std::vector<scalar> forward
    (
        const std::vector<scalar>& x
    ) const override;
};


class RescalingLayer : public Layer
{
public:

    scalar scale;
    scalar offset;

    std::vector<scalar> forward
    (
        const std::vector<scalar>& x
    ) const override;
};

std::vector<scalar> loadBin(const fileName& path, label n);

class wallModelMLP
{
private:

    std::vector<std::unique_ptr<Layer>> layers_;

public:

    wallModelMLP(const fileName& dir);

    std::vector<scalar> predict
    (
        const std::vector<scalar>& rawInputs
    ) const;
};

}

#endif