#ifndef WINDOWS_CALCENGINE_WCEUNIQUE_H
#define WINDOWS_CALCENGINE_WCEUNIQUE_H

#include <memory>

// C++ 11 does not support make_unique. Using this header file to have make_unique functionality.
// When using later C++ standards, remove this header file. Note also that some compilers might have
// make_unique in C++ 11, but that is not case for all of them. That is why we want to use this
// function

namespace fem {
	template< class T, class... Args >
	std::unique_ptr< T > make_unique( Args && ... args ) {
		return std::unique_ptr< T >( new T( std::forward< Args >( args )... ) );
	}

}

#endif //WINDOWS_CALCENGINE_WCEUNIQUE_H