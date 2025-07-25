const Title = () => {
  return (
    <header className="relative text-white">
      {/* Header Background Enhancement */}
      <div className="absolute inset-0 bg-gradient-to-b from-white/5 to-transparent" />

      <div className="relative z-10 container mx-auto px-4 py-4">
        {/* GitHub Link */}
        <div className="absolute top-4 right-4">
          <a
            href="https://github.com/BruhLemma-Yadecha/rocketmancer-web"
            target="_blank"
            rel="noopener noreferrer"
            className="inline-flex items-center glass p-2 rounded-lg hover:glass-strong transition-all duration-200 group"
            title="View source on GitHub"
            aria-label="View source code on GitHub"
          >
            <svg
              className="w-4 h-4 text-white/80 group-hover:text-white transition-colors duration-200"
              fill="currentColor"
              viewBox="0 0 20 20"
              aria-hidden="true"
            >
              <path
                fillRule="evenodd"
                d="M10 0C4.477 0 0 4.484 0 10.017c0 4.425 2.865 8.18 6.839 9.504.5.092.682-.217.682-.483 0-.237-.008-.868-.013-1.703-2.782.605-3.369-1.343-3.369-1.343-.454-1.158-1.11-1.466-1.11-1.466-.908-.62.069-.608.069-.608 1.003.07 1.531 1.032 1.531 1.032.892 1.53 2.341 1.088 2.91.832.092-.647.35-1.088.636-1.338-2.22-.253-4.555-1.113-4.555-4.951 0-1.093.39-1.988 1.029-2.688-.103-.253-.446-1.272.098-2.65 0 0 .84-.27 2.75 1.026A9.564 9.564 0 0110 4.844c.85.004 1.705.115 2.504.337 1.909-1.296 2.747-1.027 2.747-1.027.546 1.379.203 2.398.1 2.651.64.7 1.028 1.595 1.028 2.688 0 3.848-2.339 4.695-4.566 4.942.359.31.678.921.678 1.856 0 1.338-.012 2.419-.012 2.747 0 .268.18.58.688.482A10.019 10.019 0 0020 10.017C20 4.484 15.522 0 10 0z"
                clipRule="evenodd"
              />
            </svg>
          </a>
        </div>

        {/* Main Header Content */}
        <div className="text-center">
          <div className="inline-flex items-center justify-center w-12 h-12 bg-white/20 rounded-full mb-2 animate-float">
            <span className="text-2xl" role="img" aria-label="rocket">
              🚀
            </span>
          </div>
          <h1 className="text-3xl md:text-4xl font-bold bg-gradient-to-r from-white to-blue-100 bg-clip-text text-transparent animate-fade-in">
            rocketmancer
          </h1>
          <p className="text-sm md:text-base text-blue-100 mt-2 max-w-xl mx-auto animate-slide-up">
            Multi-stage rocket optimizer
          </p>
        </div>
      </div>

      {/* Decorative Background Elements */}
      <div className="absolute inset-0 overflow-hidden pointer-events-none" aria-hidden="true">
        <div className="absolute top-4 left-8 w-1 h-1 bg-white/20 rounded-full animate-pulse" />
        <div className="absolute top-8 right-12 w-0.5 h-0.5 bg-white/30 rounded-full animate-pulse delay-300" />
        <div className="absolute bottom-4 left-1/4 w-1 h-1 bg-white/25 rounded-full animate-pulse delay-700" />
        <div className="absolute bottom-6 right-1/3 w-0.5 h-0.5 bg-white/20 rounded-full animate-pulse delay-1000" />
      </div>
    </header>
  );
};

export default Title;
