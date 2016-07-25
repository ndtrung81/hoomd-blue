#ifndef _SHAPE_MOVES_H
#define _SHAPE_MOVES_H
#include <hoomd/hoomd_config.h>
#include <hoomd/saruprng.h>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
// #include "ShapeProxy.h"
#include "ShapeUtils.h"

namespace hpmc {

namespace detail{

template < typename Shape >
struct shape_param_to_vector
    {
    void operator ()(const typename Shape::param_type&, std::vector<Scalar>&) { }
    };

template <unsigned int max_verts>
struct shape_param_to_vector< ShapeConvexPolyhedron<max_verts> >
    {
    void operator ()(const typename ShapeConvexPolyhedron<max_verts>::param_type& shape, std::vector<Scalar>& params)
        {
        params.resize(shape.N*3, 0.0);
        for(size_t i = 0; i < shape.N; i++)
            {
            params[i*3]   = shape.x[i];
            params[i*3+1] = shape.y[i];
            params[i*3+2] = shape.z[i];
            }
        }
    };


// void unpack_ignore_flag(unsigned int flag, bool& stats, bool& ovrlps);

}

template < typename Shape, typename RNG>
class shape_move_function
{
public:
    shape_move_function(unsigned int ntypes) : m_determinantInertiaTensor(0), m_step_size(ntypes) {}
    shape_move_function(const shape_move_function& src) : m_determinantInertiaTensor(src.getDeterminantInertiaTensor()), m_step_size(src.getStepSize()) {}

    virtual void prepare(unsigned int timestep) = 0;
    //! construct is called at the beginning of every update()
    virtual void construct(const unsigned int&, const unsigned int&, typename Shape::param_type&, RNG&) = 0;

    //! advance whenever the proposed move is accepted.
    // virtual void advance(const unsigned int) = 0;

    //! retreat whenever the proposed move is rejected.
    virtual void retreat(const unsigned int) = 0;

    virtual Scalar getParam(size_t i ) { return 0.0; }

    virtual size_t getNumParam() { return 0; }

    Scalar getDeterminant() const { return m_determinantInertiaTensor; }

    Scalar getStepSize(const unsigned int& type_id) const { return m_step_size[type_id]; }

    void setStepSize(const unsigned int& type_id, const Scalar& stepsize) { m_step_size[type_id] = stepsize; }

protected:
    unsigned int                    m_typeid;                       // input. last typeid that was moved.
    Scalar                          m_determinantInertiaTensor;     // output
    std::vector<Scalar>             m_step_size;                    // maximum stepsize. input/output
};

template < typename Shape, typename RNG >
class python_callback_parameter_shape_move : public shape_move_function<Shape, RNG>
{
    using shape_move_function<Shape, RNG>::m_determinantInertiaTensor;
    using shape_move_function<Shape, RNG>::m_step_size;
    using shape_move_function<Shape, RNG>::m_typeid;
public:
    python_callback_parameter_shape_move(   unsigned int ntypes,
                                            boost::python::object python_function,
                                            const std::vector< std::vector<Scalar> >& params,
                                            const std::vector<Scalar>& stepsize,
                                            Scalar mixratio,
                                            bool normalized
                                        )
        :  shape_move_function<Shape, RNG>(ntypes), m_params(params), m_python_callback(python_function), m_normalized(normalized)
        {
        if(m_step_size.size() != stepsize.size())
            throw std::runtime_error("must provide a stepsize for each type");

        m_step_size = stepsize;
        m_select_ratio = fmin(mixratio, 1.0)*65535;
        boost::python::object tup = m_python_callback(m_params);
        m_determinantInertiaTensor = boost::python::extract< Scalar >(tup[1]);
        }

    void prepare(unsigned int timestep)
        {
        m_params_backup = m_params;
        m_step_size_backup = m_step_size;
        }

    void construct(const unsigned int& timestep, const unsigned int& type_id, typename Shape::param_type& shape, RNG& rng)
        {
        // gonna make the move.
        // Saru rng(m_select_ratio, m_seed, timestep);
        m_typeid = type_id;

        if(m_normalized)
            {
            for(size_t i = 0; i < m_params.size(); i++)
                {
                Scalar x = ((rng.u32() & 0xffff) < m_select_ratio) ? rng.s(fmax(-m_step_size[m_typeid], -(m_params[m_typeid][i])), fmin(m_step_size[m_typeid], (1.0-m_params[m_typeid][i]))) : 0.0;
                m_params[m_typeid][i] += x;
                }
            }
        else
            {
            // could gut this part becuase of the other functions.
            for(size_t i = 0; i < m_params[m_typeid].size() && (m_params[m_typeid].size()%3 == 0); i+=3)
                {
                if( (rng.u32()& 0xffff) < m_select_ratio )
                    {
                    Scalar x = rng.s(-1.0, 1.0);
                    Scalar y = rng.s(-1.0, 1.0);
                    Scalar z = rng.s(-1.0, 1.0);
                    Scalar mag = rng.s(0.0, m_step_size[m_typeid])/sqrt(x*x + y*y + z*z);
                    m_params[m_typeid][i] += x*mag;
                    m_params[m_typeid][i+1] += y*mag;
                    m_params[m_typeid][i+2] += z*mag;
                    }
                }
            }

        boost::python::object tup = m_python_callback(m_params[m_typeid]);
        shape = boost::python::extract< typename Shape::param_type >(tup[0]);
        m_determinantInertiaTensor = boost::python::extract< Scalar >(tup[1]);
        m_scale = Scalar(1.0);
        if(!m_normalized)
            {
            m_scale = boost::python::extract< Scalar >(tup[2]); // only extract if we have to.
            detail::shape_param_to_vector<Shape> converter;
            converter(shape, m_params[m_typeid]);
            }
        m_step_size[m_typeid] *= m_scale; // only need to scale if the parameters are not normalized
        }

    // void advance(unsigned int timestep)
    //     {
    //     // move has been accepted. nothing to do.
    //     }

    void retreat(unsigned int timestep)
        {
        // move has been rejected.
        std::swap(m_params, m_params_backup);
        std::swap(m_step_size, m_step_size_backup);
        }

    Scalar getParam(size_t k)
        {
        size_t n = 0;
        for (size_t i = 0; i < m_params.size(); i++)
            {
            size_t next = n + m_params[i].size();
            if(k < next)
                return m_params[i][k - n];
            n = next;
            }

        return 0.0; // out of range.
        }

    size_t getNumParam()
        {
        size_t n = 0;
        for (size_t i = 0; i < m_params.size(); i++)
            n += m_params[i].size();
        return n;
        }

private:
    std::vector<Scalar>                     m_step_size_backup;
    unsigned int                            m_select_ratio;     // fraction of parameters to change in each move. internal use
    Scalar                                  m_scale;            // the scale needed to keep the particle at constant volume. internal use
    std::vector< std::vector<Scalar> >      m_params_backup;    // all params are from 0,1
    std::vector< std::vector<Scalar> >      m_params;           // all params are from 0,1
    boost::python::object                   m_python_callback;  // callback that takes m_params as an argiment and returns (shape, det(I))
    bool                                    m_normalized;       // if true all parameters are restricted to (0,1)
};

template< typename Shape, typename RNG >
class constant_shape_move : public shape_move_function<Shape, RNG>
{
    using shape_move_function<Shape, RNG>::m_determinantInertiaTensor;
    using shape_move_function<Shape, RNG>::m_typeid;
public:
    constant_shape_move(const typename Shape::param_type& shape_move, Scalar detI_move) : shape_move_function<Shape, RNG>(1), m_shapeMove(shape_move), m_determinantInertiaTensorMove(detI_move)
        {
        }

    void prepare(unsigned int timestep) {}

    void construct(const unsigned int& timestep, const unsigned int& type_id, typename Shape::param_type& shape, RNG& rng)
        {
        shape = m_shapeMove;
        m_determinantInertiaTensor = m_determinantInertiaTensorMove;
        }

    // void advance(unsigned int timestep)
    //     {
    //     // move has been accepted.
    //     }

    void retreat(unsigned int timestep)
        {
        // move has been rejected.
        }

private:
    typename Shape::param_type      m_shapeMove;
    Scalar                          m_determinantInertiaTensorMove;
};

template < typename ShapeConvexPolyhedronType, typename RNG >
class convex_polyhedron_generalized_shape_move : public shape_move_function<ShapeConvexPolyhedronType, RNG>
{
    using shape_move_function<ShapeConvexPolyhedronType, RNG>::m_determinantInertiaTensor;
    using shape_move_function<ShapeConvexPolyhedronType, RNG>::m_step_size;
    using shape_move_function<ShapeConvexPolyhedronType, RNG>::m_typeid;
public:
    convex_polyhedron_generalized_shape_move(
                                            unsigned int ntypes,
                                            Scalar stepsize,
                                            Scalar mixratio,
                                            Scalar volume
                                        ) : shape_move_function<ShapeConvexPolyhedronType, RNG>(ntypes), m_volume(volume)
        {
                    // if(m_step_size.size() != stepsize.size())
                    //     throw std::runtime_error("must provide a stepsize for each type");

        m_determinantInertiaTensor = 1.0;
        m_scale = 1.0;
        m_step_size.clear();
        m_step_size.resize(ntypes, stepsize);
        m_select_ratio = fmin(mixratio, 1.0)*65535;
        }

    void prepare(unsigned int timestep)
        {
        m_step_size_backup = m_step_size;
        }

    void construct(const unsigned int& timestep, const unsigned int& type_id, typename ShapeConvexPolyhedronType::param_type& shape, RNG& rng)
        {
        m_typeid = type_id;
        // mix the shape.
        for(size_t i = 0; i < shape.N; i++)
            {
            if( (rng.u32()& 0xffff) < m_select_ratio )
                {
                Scalar x = rng.s(-1.0, 1.0);
                Scalar y = rng.s(-1.0, 1.0);
                Scalar z = rng.s(-1.0, 1.0);
                Scalar mag = rng.s(0.0, m_step_size[m_typeid])/sqrt(x*x + y*y + z*z);
                shape.x[i] += x*mag;
                shape.y[i] += y*mag;
                shape.z[i] += z*mag;
                }
            }

        detail::ConvexHull convex_hull(shape); // compute the convex_hull.
        convex_hull.compute();
        detail::mass_properties<ShapeConvexPolyhedronType> mp(convex_hull.getPoints(), convex_hull.getFaces());
        Scalar volume = mp.getVolume();
        vec3<Scalar> centroid = mp.getCenterOfMass();
        m_scale = pow(m_volume/volume, 1.0/3.0);
        Scalar rsq = 0.0;
        std::vector< vec3<Scalar> > points(shape.N);
        for(size_t i = 0; i < shape.N; i++)
            {
            shape.x[i] -= centroid.x;
            shape.x[i] *= m_scale;
            shape.y[i] -= centroid.y;
            shape.y[i] *= m_scale;
            shape.z[i] -= centroid.z;
            shape.z[i] *= m_scale;
            vec3<Scalar> vert(shape.x[i], shape.y[i], shape.z[i]);
            rsq = fmax(rsq, dot(vert, vert));
            points[i] = vert;
            }
        detail::mass_properties<ShapeConvexPolyhedronType> mp2(points, convex_hull.getFaces());
        m_determinantInertiaTensor = mp2.getDeterminant();
        shape.diameter = 2.0*sqrt(rsq);
        m_step_size[m_typeid] *= m_scale; // only need to scale if the parameters are not normalized
        }

    // void advance(unsigned int timestep)
    //     {
    //     // nothing to do.
    //     }

    void retreat(unsigned int timestep)
        {
        // move has been rejected.
        std::swap(m_step_size, m_step_size_backup);
        }

private:
    std::vector<Scalar>     m_step_size_backup;
    unsigned int            m_select_ratio;
    Scalar                  m_scale;
    Scalar                  m_volume;
};

template <class Shape, class RNG>
struct shear
    {
    shear(Scalar) {}
    void operator() (typename Shape::param_type& param, RNG& rng)
        {
        throw std::runtime_error("shear is not implemented for this shape.");
        }
    };

template <class Shape, class RNG>
struct scale
    {
    bool isotropic;
    scale(bool iso = true) : isotropic(iso) {}
    void operator() (typename Shape::param_type& param, RNG& rng)
        {
        throw std::runtime_error("scale is not implemented for this shape.");
        }
    };


template <unsigned int max_verts, class RNG>
struct shear< ShapeConvexPolyhedron<max_verts>, RNG >
    {
    Scalar shear_max;
    shear(Scalar smax) : shear_max(smax) {}
    void operator() (typename ShapeConvexPolyhedron<max_verts>::param_type& param, RNG& rng)
        {
        Scalar gamma = rng.s(-shear_max, shear_max), gammaxy = 0.0, gammaxz = 0.0, gammayz = 0.0, gammayx = 0.0, gammazx = 0.0, gammazy = 0.0;
        int dim = int(6*rng.s(0.0, 1.0));
        if(dim == 0) gammaxy = gamma;
        else if(dim == 1) gammaxz = gamma;
        else if(dim == 2) gammayz = gamma;
        else if(dim == 3) gammayx = gamma;
        else if(dim == 4) gammazx = gamma;
        else if(dim == 5) gammazy = gamma;
        Scalar dsq = 0.0;
        for(unsigned int i = 0; i < param.N; i++)
            {
            param.x[i] = param.x[i] + param.y[i]*gammaxy + param.z[i]*gammaxz;
            param.y[i] = param.x[i]*gammayx + param.y[i] + param.z[i]*gammayz;
            param.z[i] = param.x[i]*gammazx + param.y[i]*gammazy + param.z[i];
            vec3<Scalar> vert( param.x[i], param.y[i], param.z[i]);
            dsq = fmax(dsq, dot(vert, vert));
            }
        param.diameter = 2.0*sqrt(dsq);
        }
    };

template <unsigned int max_verts, class RNG>
struct scale< ShapeConvexPolyhedron<max_verts>, RNG >
    {
    bool isotropic;
    Scalar scale_min;
    Scalar scale_max;
    scale(Scalar movesize, bool iso = true) : isotropic(iso)
    {
    if(movesize < 0.0 || movesize > 1.0)
        {
        movesize = 0.0;
        }
    scale_max = (1.0+movesize);
    scale_min = 1.0/scale_max;
    }
    void operator() (typename ShapeConvexPolyhedron<max_verts>::param_type& param, RNG& rng)
        {
        Scalar sx, sy, sz;
        Scalar s = rng.s(scale_min, scale_max);
        sx = sy = sz = s;
        if(!isotropic)
            {
            sx = sy = sz = 1.0;
            Scalar dim = rng.s(0.0, 1.0);
            if (dim < 1.0/3.0) sx = s;
            else if (dim < 2.0/3.0) sy = s;
            else sz = s;
            }
        for(unsigned int i = 0; i < param.N; i++)
            {
            param.x[i] *= sx;
            param.y[i] *= sy;
            param.z[i] *= sz;
            }
        param.diameter *= s;
        }
    };

template < class RNG >
class scale< ShapeEllipsoid, RNG >
{
    const Scalar m_v;
    const Scalar m_v1;
    const Scalar m_min;
    const Scalar m_max;
public:
    scale(Scalar movesize, bool) : m_v(1.0), m_v1(M_PI*4.0/3.0), m_min(-movesize), m_max(movesize) {}
    void operator ()(ShapeEllipsoid::param_type& param, RNG& rng)
        {
        Scalar lnx = log(param.x/param.y);
        Scalar dx = rng.s(m_min, m_max);
        Scalar x = fast::exp(lnx+dx);
        Scalar b = pow(m_v/m_v1/x, 1.0/3.0);

        param.x = x*b;
        param.y = b;
        param.z = b;
        }
};



template<class Shape, class RNG>
class elastic_shape_move_function : public shape_move_function<Shape, RNG>
{
    // using shape_move_function<Shape, RNG>::m_shape;
    using shape_move_function<Shape, RNG>::m_determinantInertiaTensor;
    using shape_move_function<Shape, RNG>::m_step_size;
    using shape_move_function<Shape, RNG>::m_typeid;
    // using shape_move_function<Shape, RNG>::m_scale;
    // using shape_move_function<Shape, RNG>::m_select_ratio;
public:
    elastic_shape_move_function(
                                    unsigned int ntypes,
                                    const Scalar& stepsize,
                                    Scalar move_ratio
                                ) : shape_move_function<Shape, RNG>(ntypes)
        {
        m_select_ratio = fmin(move_ratio, 1.0)*65535;
        m_step_size.resize(ntypes, stepsize);
        std::fill(m_step_size.begin(), m_step_size.end(), stepsize);
        }

    void prepare(unsigned int timestep) { /* Nothing to do. */ }
    //! construct is called at the beginning of every update()
    void construct(const unsigned int& timestep, const unsigned int& type_id, typename Shape::param_type& shape, RNG& rng)
        {
        unsigned int move_type_select = rng.u32() & 0xffff;
        // Saru rng(m_select_ratio, m_seed, timestep);
        if( move_type_select < m_select_ratio)
            {
            scale<Shape, RNG> move(m_step_size[type_id], false);
            move(shape, rng); // always make the move
            }
        else
            {
            shear<Shape, RNG> move(m_step_size[type_id]);
            move(shape, rng); // always make the move
            }
        detail::mass_properties<Shape> mp(shape); // this could be slow for some shapes.
        m_determinantInertiaTensor = mp.getDeterminant();
        }

    //! advance whenever the proposed move is accepted.
    // void advance(unsigned int timestep){ /* Nothing to do. */ }

    //! retreat whenever the proposed move is rejected.
    void retreat(unsigned int timestep){ /* Nothing to do. */ }

protected:
    unsigned int            m_select_ratio;
};

template<class Shape>
class ShapeLogBoltzmannFunction
{
public:
    virtual Scalar operator()(const unsigned int& N, const typename Shape::param_type& shape_new, const Scalar& inew, const typename Shape::param_type& shape_old, const Scalar& iold) { throw std::runtime_error("not implemented"); return 0.0;}
};

template<class Shape>
class AlchemyLogBoltzmannFunction : public ShapeLogBoltzmannFunction<Shape>
{
public:
    virtual Scalar operator()(const unsigned int& N, const typename Shape::param_type& shape_new, const Scalar& inew, const typename Shape::param_type& shape_old, const Scalar& iold)
        {
        return (Scalar(N)/Scalar(2.0))*log(inew/iold);
        }
};

class EllipsoidSpringPRL : public ShapeLogBoltzmannFunction<ShapeEllipsoid>
{
    Scalar m_k;
public:
    EllipsoidSpringPRL(Scalar k) : m_k(k) {}
    Scalar operator()(const unsigned int& N, const ShapeEllipsoid::param_type& shape_new, const Scalar& inew, const ShapeEllipsoid::param_type& shape_old, const Scalar& iold)
        {
        Scalar x_new = shape_new.x/shape_new.y;
        Scalar x_old = shape_old.x/shape_old.y;
        return m_k*(log(x_old)*log(x_old) - log(x_new)*log(x_new)); // -\beta dH
        }
};

template<class ConvexPolyhedronShape>
class ConvexPolyhedronSpring : public ShapeLogBoltzmannFunction<ConvexPolyhedronShape>
{
    Scalar m_k;
    typename ConvexPolyhedronShape::param_type m_reference_shape;
public:
    ConvexPolyhedronSpring(Scalar k, const typename ConvexPolyhedronShape::param_type& ref ) : m_k(k), m_reference_shape(ref) {}
    Scalar operator()(const unsigned int& N, const typename ConvexPolyhedronShape::param_type& shape_new, const Scalar& inew, const typename ConvexPolyhedronShape::param_type& shape_old, const Scalar& iold)
        {
        AlchemyLogBoltzmannFunction<ConvexPolyhedronShape> fn;
        Scalar U_old = 0.0, U_new = 0.0;
        for(unsigned int i = 0; i < m_reference_shape.N; i++)
            {
            vec3<Scalar> v_old(shape_old.x[i], shape_old.y[i], shape_old.z[i]), v_new(shape_new.x[i], shape_new.y[i], shape_new.z[i]), v_ref(m_reference_shape.x[i], m_reference_shape.y[i], m_reference_shape.z[i]), dro, drn;
            dro = v_old - v_ref;
            drn = v_new - v_ref;
            U_old += m_k*dot(dro, dro);
            U_new += m_k*dot(drn, drn);
            }
        return (U_old - U_new) + fn(N, shape_new, inew, shape_old, iold); // -\beta dH
        }
};

//** Python export functions and additional classes to wrap the move and boltzmann interface.
//**
//**
//**
//**
//! Wrapper class for wrapping pure virtual methods
template<class Shape, class RNG>
class shape_move_function_wrap : public shape_move_function<Shape, RNG>, public boost::python::wrapper< shape_move_function<Shape, RNG> >
    {
    public:
        //! Constructor
        shape_move_function_wrap(unsigned int ntypes) : shape_move_function<Shape, RNG>(ntypes) {}
        void prepare(unsigned int timestep)
            {
            this->get_override("prepare")(timestep);
            }

        void construct(const unsigned int& timestep, const unsigned int& type_id, typename Shape::param_type& shape, RNG& rng)
            {
            this->get_override("construct")(timestep, type_id, shape, rng);
            }

        // void advance(unsigned int timestep)
        //     {
        //     this->get_override("advance")(timestep);
        //     }

        void retreat(unsigned int timestep)
            {
            this->get_override("retreat")(timestep);
            }
    };

template<class Shape>
void export_ShapeMoveInterface(std::string name);

template<class Shape>
void export_ElasticShapeMove(std::string name);

template< typename Shape >
void export_ShapeLogBoltzmann(const std::string& name);

void export_LogBoltzmannEllipsoidSpring(const std::string& name);

template<class Shape>
void export_LogBoltzmannConvexPolyhedronSpring(const std::string& name);

template<class Shape>
void export_AlchemyLogBoltzmannFunction(const std::string& name);

template<class Shape>
void export_ConvexPolyhedronGeneralizedShapeMove(const std::string& name);

}

#endif
