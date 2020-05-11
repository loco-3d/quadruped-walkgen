
#ifndef __mpc_controller_macros_hpp__
#define __mpc_controller_macros_hpp__

#define MPCC_BOOST_MAKE_SHARED_PTR(classtype, ptrname, ...) boost::shared_ptr<classtype> ptrname = boost::make_shared<classtype>(__VA_ARGS__)



#endif //  __mpc_controller_macros_hpp__
