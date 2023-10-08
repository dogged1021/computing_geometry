
namespace ComputingGeometry
{
    template <typename T>
    class GeometryBase
    {
    public:
        virtual GeometryBase() = 0;
        virtual ~GeometryBase() = 0;

    private:
        T size;
    };

} // namespace ComputingGeometry
