// Allow expression of type lambda * V by forwarding to V * lambda
template <typename Derived, typename TVal, typename TIdx>
BinaryOperation<operation::type::scalar_product, Derived, TVal>
operator* (const TVal& alpha, const DMatrixBase<Derived, TVal, TIdx>& base)
{
    return BinaryOperation<operation::type::scalar_product,
           Derived,
           TVal>(base.derived(), alpha);
}
