namespace Zee {

template <typename TVal, typename TIdx>
class DStreamingVector;

template <typename TVal, typename TIdx>
class DStreamingSparseMatrix :
    public DSparseMatrix<TVal, TIdx>
{
    public:
        DStreamingSparseMatrix(
                std::shared_ptr<UnpainBase::Center<TIdx>> center,
                std::string file,
                TIdx procs = 0) :
            DSparseMatrix<TVal, TIdx>(center, file, procs)
        { }

        BinaryOperation<operation::type::product,
            DStreamingSparseMatrix<TVal, TIdx>,
            DStreamingVector<TVal, TIdx>>
            operator* (const DStreamingVector<TVal, TIdx>& rhs) const
        {
            return BinaryOperation<operation::type::product,
                DStreamingSparseMatrix<TVal, TIdx>,
                DStreamingVector<TVal, TIdx>>(*this, rhs);
        }
};

template <typename TVal, typename TIdx>
class DStreamingVector :
    public DVector<TVal, TIdx>
{
    public:
        DStreamingVector(std::shared_ptr<UnpainBase::Center<TIdx>> center, TIdx n, TVal defaultValue = (TVal)0)
            : DVector<TVal, TIdx>(center, n, defaultValue)
        {
        }

        using DVector<TVal, TIdx>::operator=;

        void operator= (const BinaryOperation<operation::type::product,
            DStreamingSparseMatrix<TVal, TIdx>,
            DStreamingVector<TVal, TIdx>>& op)
        {
            // TODO
            ZeeLogInfo << "Epiphany SpMV" << endLog;
        }
};

} // namespace Zee
