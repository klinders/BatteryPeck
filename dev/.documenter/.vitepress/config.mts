import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import { mathjaxPlugin } from './mathjax-plugin'
import { juliaReplTransformer } from './julia-repl-transformer'
import footnote from "markdown-it-footnote";
import path from 'path'

const mathjax = mathjaxPlugin()

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: '/BatteryToolkit/',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: [
{ text: 'Home', link: '/index' },
{ text: 'Getting Started', collapsed: false, items: [
{ text: 'Quick Start', link: '/guide/quickstart' },
{ text: 'Parameter Sets', link: '/guide/parameters' },
{ text: 'Experiments', link: '/guide/experiments' },
{ text: 'Examples', link: '/guide/examples' }]
 },
{ text: 'Models', collapsed: false, items: [
{ text: 'SPMe Cell Model', link: '/models/spme' },
{ text: 'Pack Models', link: '/models/pack-models' },
{ text: 'Side Reactions', link: '/models/side-reactions' }]
 },
{ text: 'API Reference', collapsed: false, items: [
{ text: 'Parameters', link: '/api/parameters' },
{ text: 'Cell Models', link: '/api/cellmodels' },
{ text: 'Pack Models', link: '/api/packmodels' },
{ text: 'Experiments', link: '/api/experiments' },
{ text: 'Solvers', link: '/api/solvers' }]
 },
{ text: 'Advanced', collapsed: false, items: [
{ text: 'Finite Volume Method', link: '/Finite Volume Method/index' }]
 },
{ text: 'References', link: '/references' }
]
,
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/BatteryToolkit/',// TODO: replace this in makedocs!
  title: 'BatteryToolkit.jl',
  description: 'Documentation for BatteryToolkit',
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../1', // This is required for MarkdownVitepress to work correctly...
  head: [
    
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  
  markdown: {
    codeTransformers: [juliaReplTransformer()],
    config(md) {
      md.use(tabsMarkdownPlugin);
      md.use(footnote);
      mathjax.markdownConfig(md);
    },
    theme: {
      light: "github-light",
      dark: "github-dark"
    },
  },
  vite: {
    plugins: [
      mathjax.vitePlugin,
    ],
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('/BatteryToolkit'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    optimizeDeps: {
      exclude: [ 
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ], 
    }, 
    ssr: { 
      noExternal: [ 
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ], 
    },
  },
  themeConfig: {
    outline: 'deep',
    logo: { src: '/logo.svg', width: 24, height: 24},
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
{ text: 'Home', link: '/index' },
{ text: 'Getting Started', collapsed: false, items: [
{ text: 'Quick Start', link: '/guide/quickstart' },
{ text: 'Parameter Sets', link: '/guide/parameters' },
{ text: 'Experiments', link: '/guide/experiments' },
{ text: 'Examples', link: '/guide/examples' }]
 },
{ text: 'Models', collapsed: false, items: [
{ text: 'SPMe Cell Model', link: '/models/spme' },
{ text: 'Pack Models', link: '/models/pack-models' },
{ text: 'Side Reactions', link: '/models/side-reactions' }]
 },
{ text: 'API Reference', collapsed: false, items: [
{ text: 'Parameters', link: '/api/parameters' },
{ text: 'Cell Models', link: '/api/cellmodels' },
{ text: 'Pack Models', link: '/api/packmodels' },
{ text: 'Experiments', link: '/api/experiments' },
{ text: 'Solvers', link: '/api/solvers' }]
 },
{ text: 'Advanced', collapsed: false, items: [
{ text: 'Finite Volume Method', link: '/Finite Volume Method/index' }]
 },
{ text: 'References', link: '/references' }
]
,
    sidebarDrawer: false,
    editLink: { pattern: "https://github.com/klinders/BatteryToolkit/edit/main/docs/src/:path" },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/klinders/BatteryToolkit' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
