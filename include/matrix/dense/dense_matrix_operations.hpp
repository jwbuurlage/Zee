template <typename TVal, typename TIdx>
DVector<TVal, TIdx> perform_operation(
        BinaryOperation<operation::type::product,
        DSparseMatrix<TVal, TIdx>,
        DVector<TVal, TIdx>> op)
{
    using TImage = typename DSparseMatrix<TVal, TIdx>::image_type;

    const auto& A = op.getLHS();
    const auto& v = op.getRHS();
    DVector<TVal, TIdx> u{A.getRows(), 0.0};
    const auto p = A.getProcs();

    ZeeAssert(A.getCols() == v.size());

    std::mutex writeMutex;
    Barrier<TIdx> barrier(p);

    A.compute([&] (std::shared_ptr<TImage> submatrixPtr,
                TIdx s) {
        std::map<TIdx, TVal> u_s;
        for (const auto& triplet : *submatrixPtr) {
             u_s[triplet.row()] += triplet.value() * v[triplet.col()];
        }

        // if RHS is a reference to *this then we can only reset here,
        // so we use a barrier sync
        barrier.sync();

        for (TIdx i = s; i < u.size(); i += p) {
            u[i] = 0;
        }

        barrier.sync();

        for (const auto& pKeyValue : u_s) {
            std::lock_guard<std::mutex> lock(writeMutex);
            u[pKeyValue.first] += pKeyValue.second;
        }
    });

    return u;
}
