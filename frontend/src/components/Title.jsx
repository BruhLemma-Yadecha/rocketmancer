const Title = () => {
  return (
    <header className="relative bg-gradient-to-r from-blue-600 via-purple-600 to-indigo-600 text-white">
      <div className="absolute inset-0 bg-black/10"></div>
      <div className="relative z-10 container mx-auto px-4 py-8">
        <div className="text-center">
          <div className="inline-flex items-center justify-center w-16 h-16 bg-white/20 rounded-full mb-4 animate-float">
            <span className="text-3xl">🚀</span>
          </div>
          <h1 className="text-4xl md:text-5xl lg:text-6xl font-bold bg-gradient-to-r from-white to-blue-100 bg-clip-text text-transparent animate-fade-in">
            rocketmancer
          </h1>
          <p className="text-lg md:text-xl text-blue-100 mt-4 max-w-2xl mx-auto animate-slide-up">
            Advanced rocket optimization and design platform
          </p>
          <div className="flex items-center justify-center mt-6 space-x-2">
            <div className="h-1 w-8 bg-white/60 rounded-full"></div>
            <div className="h-1 w-4 bg-white/40 rounded-full"></div>
            <div className="h-1 w-2 bg-white/30 rounded-full"></div>
          </div>
        </div>
      </div>
      
      {/* Decorative elements */}
      <div className="absolute top-0 left-0 w-full h-full overflow-hidden">
        <div className="absolute top-10 left-10 w-2 h-2 bg-white/20 rounded-full animate-pulse"></div>
        <div className="absolute top-20 right-20 w-1 h-1 bg-white/30 rounded-full animate-pulse delay-300"></div>
        <div className="absolute bottom-10 left-1/4 w-1.5 h-1.5 bg-white/25 rounded-full animate-pulse delay-700"></div>
        <div className="absolute bottom-20 right-1/3 w-1 h-1 bg-white/20 rounded-full animate-pulse delay-1000"></div>
      </div>
    </header>
  );
};

export default Title;
