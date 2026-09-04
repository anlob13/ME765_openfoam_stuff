#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>

#include "scalar.H"
#include "fileName.H"
#include "error.H"

#include "mlForwardPass_full.H"

namespace Foam
{

// --- Activation functions ---
inline scalar relu(scalar x) { return std::max<scalar>(0, x); }
inline scalar tanhAct(scalar x) { return std::tanh(x); }
inline scalar linear(scalar x) { return x; }

std::vector<scalar> DenseLayer::forward
(
    const std::vector<scalar>& x
) const
{
    if (label(x.size()) != inDim)
    {
        FatalErrorInFunction
            << "Input size " << x.size()
            << " does not match expected " << inDim
            << exit(FatalError);
    }

    std::vector<scalar> y(outDim, 0.0);

    for (label j=0; j<outDim; ++j)
    {
        scalar sum = b[j];

        for (label i=0; i<inDim; ++i)
        {
            sum += x[i]*W[i*outDim + j];
        }

        y[j] = activation(sum);
    }

    return y;
}


std::vector<scalar> NormalizationLayer::forward
(
    const std::vector<scalar>& x
) const
{
    if (x.size() != mean.size())
    {
        FatalErrorInFunction
            << "Normalization dimension mismatch."
            << exit(FatalError);
    }

    std::vector<scalar> y(x.size());

    for (label i=0; i<label(x.size()); ++i)
    {
        y[i] =
            (x[i] - mean[i])/
            std::sqrt(variance[i]);
    }

    return y;
}


std::vector<scalar> RescalingLayer::forward
(
    const std::vector<scalar>& x
) const
{
    std::vector<scalar> y(x.size());

    for (label i=0; i<label(x.size()); ++i)
    {
        y[i] = x[i]*scale + offset;
    }

    return y;
}

// --- Load a flat binary blob into a vector ---
std::vector<scalar> loadBin(const fileName& path, label n)
{
    std::vector<scalar> data(n);
    std::ifstream f(path, std::ios::binary);
    if (!f) FatalErrorInFunction << "Cannot open " << path << exit(FatalError);
    f.read(reinterpret_cast<char*>(data.data()), n * sizeof(scalar));
    return data;
}

wallModelMLP::wallModelMLP(const fileName& dir)
:
    layers_()
{
    std::ifstream arch(dir/"layer_metadata.txt");

    if (!arch)
    {
        FatalErrorInFunction
            << "Cannot open metadata file "
            << dir/"layer_metadata.txt"
            << exit(FatalError);
    }

    int layer = 0;

    std::string type;

    while (arch >> type)
    {
        if (type == "Input")
        {
            label size;
            arch >> size;
            layer++;
        }

        else if (type == "Normalization")
        {
            label size;
            arch >> size;

            auto norm = std::make_unique<NormalizationLayer>();

            norm->mean =
                loadBin(dir/"layer"+Foam::name(layer)+"_mean.bin", size);

            norm->variance =
                loadBin(dir/"layer"+Foam::name(layer)+"_var.bin", size);

            layers_.push_back(std::move(norm));

            layer++;
        }

        else if (type == "Dense")
        {
            int inSize, outSize;
            std::string actName;

            arch >> inSize >> outSize >> actName;

            scalar (*act)(scalar) = nullptr;

            if (actName == "tanh")
                act = tanhAct;
            else if (actName == "relu")
                act = relu;
            else if (actName == "linear")
                act = linear;
            else
                FatalError
                    << "Unknown activation "
                    << actName
                    << exit(FatalError);


            layers_.push_back
            (
                std::make_unique<DenseLayer>
                (
                    inSize,
                    outSize,
                    loadBin
                    (
                        dir/("layer"+Foam::name(layer)+"_W.bin"),
                        inSize*outSize
                    ),
                    loadBin
                    (
                        dir/("layer"+Foam::name(layer)+"_b.bin"),
                        outSize
                    ),
                    act
                )
            );

            layer++;
        }

        else if (type == "Rescaling")
        {
            auto scale =
                std::make_unique<RescalingLayer>();

            std::ifstream f(dir/"layer"+Foam::name(layer)+"_rescaling.txt");

            f >> scale->scale;
            f >> scale->offset;

            layers_.push_back(std::move(scale));

            layer++;
        }
    }

}


std::vector<scalar> wallModelMLP::predict
(
    const std::vector<scalar>& rawInputs
) const
{
    std::vector<scalar> x(rawInputs.size());
    for (size_t i = 0; i < x.size(); ++i)
        x[i] = (rawInputs[i]); // - inMean_[i]) / inStd_[i];

    for (const auto& layer : layers_)
    {
        x = layer->forward(x);
    }

    return x; 
}


}