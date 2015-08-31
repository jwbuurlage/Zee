namespace Zee {

template <typename TVal, typename TIdx>
class DStreamingSparseMatrix :
    public DSparseMatrixBase<
        DStreamingSparseMatrix<TVal, TIdx>, TVal, TIdx>
{
    public:
        using Base = DSparseMatrixBase<
            DStreamingSparseMatrix<TVal, TIdx>, TVal, TIdx>;

        DStreamingSparseMatrix(std::string file,
                TIdx procs = 0) :
            Base(file, procs)
        { }
};

template <typename TVal = default_scalar_type,
         typename TIdx = default_index_type>
class DStreamingVector
    : public DVectorBase<DStreamingVector<TVal, TIdx>, TVal, TIdx>
{
    public:
        using Base = DVectorBase<DStreamingVector<TVal, TIdx>, TVal, TIdx>;
        using Base::operator=;

        DStreamingVector(TIdx size, TVal defaultValue = 0.0)
            : Base(size, defaultValue)
        { }
};

template <typename TVal, typename TIdx>
DStreamingVector<TVal, TIdx> perform_operation(
        BinaryOperation<operation::type::product,
        DStreamingSparseMatrix<TVal, TIdx>,
        DStreamingVector<TVal, TIdx>> op)
{
    const auto& A = op.getLHS();
    const auto& x = op.getRHS();
    DStreamingVector<TVal, TIdx> y(A.getRows(), 1.0);

    ZeeLogInfo << "SpMV on Epiphany" << endLog;
    ZeeLogVar(A.nonZeros());
    ZeeLogVar(x.size());
    return y;
}

} // namespace Zee
